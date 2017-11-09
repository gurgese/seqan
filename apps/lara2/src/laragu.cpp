// ==========================================================================
//               LaRAgu - Lagrangian Relaxation Aligner GU
// ==========================================================================
// Copyright (c) 2015-2016, Gianvito Urgese
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Gianvito Urgese nor the names of its contributors
//       may be used to endorse or promote products derived from this software
//       without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL GIANVITO URGESE OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Gianvito Urgese <gianvito.urgese@polito.it>
// ==========================================================================
// This file contains the seqan_laragu application.
// ==========================================================================

#define SEQAN_LARAGU

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <limits>
#include <sstream>
//#include <omp.h>
#include <ctime>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/arg_parse.h>
#include <seqan/store.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/align_rna.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>
#include <seqan/rna_io.h>


// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// defines all the constants used in the app
#include "data_types.h"
#include "option.h"
#include "lara_core.h"
#include "alignment.h"
#include "lemon_graph.h"
#include "tcoffee_interface.h"
#include "lara_io.h"

using namespace seqan;

// ----------------------------------------------------------------------------
// Function main()
// ----------------------------------------------------------------------------

int main (int argc, char const ** argv)
{
    // Check for absence of options.
    if (argc == 1)
    {
        std::cout << "Type " << argv[0] << " --help to get the parameters table (-i option is mandatory)." << std::endl;
        return 1;
    }

    // Parse options.
    ArgumentParser argumentParser;
    LaraOptions options;
    setupArgumentParser(argumentParser, options);
    ArgumentParser::ParseResult res = parse(options, argumentParser, argc, argv); // Fill the LaraOptions
    if (res != ArgumentParser::ParseResult::PARSE_OK)
        return res == ArgumentParser::ParseResult::PARSE_ERROR ? 1 : 0;

    // Read input files.
    // If one input file is given, then build unique pairs of the input sequences.
    // If two input files are given, then build the cross product of all sequences from the first file
    //     with all sequences from the second file.
    RnaStructContentsPair filecontents;
    readRnaFile(filecontents.first, options.inFile, options);
    readRnaFile(filecontents.second, options.inFileRef, options);
    _VV(options, "Read " << length(filecontents.first.records) << " and "
                         << length(filecontents.second.records) << " records from input files.");

    // If not present, compute the weighted interaction edges using ViennaRNA functions.
    computeMissingInteractions(filecontents.first.records, options);
    computeMissingInteractions(filecontents.second.records, options);
    _VVV(options, getEbpseqString(filecontents.first) << getEbpseqString(filecontents.second));


    // A StringSet of SeqAn-style alignments is created
    StringSet<RnaAlignment> alignments;
    RnaAlignmentTraitsVector aliTraitsVect;
    prepareAlignmentVector(alignments, aliTraitsVect, filecontents);
    _VV(options, "Number of pairwise alignments to be computed: " << length(alignments));

    // Bitvector that expresses whether an alignment is finished (i.e. bound difference is very small)
    std::vector<bool> isOptimal;
    isOptimal.resize(length(alignments), false);

    for (RnaAlignmentTraits & traits : aliTraitsVect)
    {
        // The lambda structs are referenced by the index of the first sequence.
        resize(traits.lamb, numVertices(traits.bppGraphH.inter));
        // Allocate capacity for mask: the length of the shorter sequence is the maximum number of expressed lines.
        reserve(traits.mask, std::min(numVertices(traits.bppGraphH.inter), numVertices(traits.bppGraphV.inter)));

        traits.structScore.score_matrix = options.laraScoreMatrix;
        traits.structScore.lamb = & traits.lamb;
        traits.stepSizeScaling = options.my;
    }

    // Create the alignment data structure that will host the alignments with small difference between bounds
    // Move all finished alignments to goldRnaAligns, such that they are not further processed
    RnaAlignmentTraitsVector goldRnaAligns;

    // Receives the alignment scores
    String<TScoreValue> resultsSimd;

    firstSimdAlignsGlobalLocal(resultsSimd, alignments, options);
    _VV(options, "\ninitial sequence alignment (score " << resultsSimd[0] << "):\n" << alignments[0]);

    for (unsigned i = 0; i < length(alignments); ++i)
    {
        RnaAlignmentTraits & traits = aliTraitsVect[i];
        createMask(traits, alignments[i], options);
        createNewLambdaLines(traits, options.useOppositLineUB, 0);
    }


    // ITERATIONS OF ALIGNMENT AND MWM
    for (unsigned iter = 0; iter < options.iterations && length(alignments) > 0; ++iter)
    {
        bool foundAnOptimalAlignment = false;
        simdAlignsGlobalLocal(resultsSimd, alignments, aliTraitsVect, options);
        _VV(options, "\nalignment in iteration " << iter << " (score " << resultsSimd[0] << "):\n" << alignments[0]);

        //#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignments); ++i)
        {
            RnaAlignmentTraits & traits = aliTraitsVect[i];
            bool changedMask = createMask(traits, alignments[i], options);

            traits.upperBound = resultsSimd[i];
            traits.bestUpperBound = std::min(traits.bestUpperBound, traits.upperBound);
//            std::cerr << "Saved Upper Bound (alignment Score)" << std::endl;

            if (changedMask || i == 0)
            {

                createNewLambdaLines(traits, options.useOppositLineUB, iter);
//                std::cerr << "Included in Lambda vector the new alignment lines" << std::endl;

                // The MWM is computed to fill the LowerBound
                if (options.lowerBoundMethod == LBLEMONMWM) {
                    TMapVect lowerBound4Lemon;
                    lowerBound4Lemon.resize(numVertices(traits.bppGraphH.inter)); //TODO check this
                    //std::cout << alignments[i] << std::endl;
                    computeLowerBound(traits, & lowerBound4Lemon);
//                computeBounds(traits, & lowerBound4Lemon); // weightLineVect receives seq indices of best pairing
//                computeUpperBoundScore(traits); // upperBound = sum of all probability lines
                    myLemon::computeLowerBoundScore(lowerBound4Lemon, traits);
                    traits.lowerBound = traits.lowerLemonBound.mwmPrimal + traits.sequenceScore;
                    // traits.slm = traits.slm - (traits.lowerLemonBound.mwmCardinality * 2);
                    _VV(options, "Computed maximum weighted matching using the LEMON library.");
                    traits.bestLowerBound = std::max(traits.bestLowerBound, traits.lowerBound);
                }
                else if (options.lowerBoundMethod == LBAPPROXMWM)
                {   // TODO Verify the procedure! Approximation of MWM is computed to fill the LowerBound
//                    computeBounds(traits, NULL);  // TODO verify this call and reimplement the approximation if better than gready
                    //TODO rewrite this function without weightLineVect
                } else if (options.lowerBoundMethod == LBLINEARTIMEMWM) // TODO Verify the procedure! using greedy algorithm
                {
                    TMapVect lowerBound4Lemon;
                    lowerBound4Lemon.resize(length(traits.mask));
//                    computeBounds(traits, &lowerBound4Lemon); //TODO rewrite this function without weightLineVect
                    computeLowerBoundGreedy(lowerBound4Lemon, traits);
                    traits.lowerBound = traits.lowerGreedyBound;
                    // traits.slm = traits.slm - (traits.lowerLemonBound.mwmCardinality * 2);
                }
            }
            _VV(options,"best l/u bounds: " << traits.bestLowerBound << " / " << traits.bestUpperBound);

            if ((traits.bestUpperBound - traits.bestLowerBound < options.epsilon))
            {
                // alignment is finished
                isOptimal[i] = true;
                foundAnOptimalAlignment = true;
                _VV(options, "Computation for alignment " << i << " stops in iteration "
                                                          << iter << " and the bestAlignMinBounds is returned.");
            }
            else
            {
                seqan::String<std::pair <unsigned, unsigned> > listUnclosedLoopMask;
                // Compute slm factor
                computeNumberOfSubgradients(traits, listUnclosedLoopMask);

                // save previous step size
                double const prev_stepSize = traits.stepSize;
                //  Compute the step size for the Lambda update
                if (traits.numberOfSubgradients > 0)
                    traits.stepSize = traits.stepSizeScaling * ((traits.bestUpperBound - traits.bestLowerBound) / traits.numberOfSubgradients);
                else
                    traits.stepSize = 0;


                //  Check the number of non decreasing iterations
                if (prev_stepSize <= traits.stepSize)
                {
                    ++traits.nonDecreasingIterations;
                    _VV(options, "nonDecreasingIterations = " << traits.nonDecreasingIterations);
                }
                else
                {
                    //TODO evaluate if the reset of this value is the right strategy wrt. the decremental solution
                    traits.nonDecreasingIterations = 0u;
                }

                // there was nothing going on in the last couple of iterations, half traits.stepSizeScaling therefore
                if (traits.nonDecreasingIterations >= options.nonDecreasingIterations)
                {
                    // this value in case of decreasing stepsize
                    // (an opposite mechanism or a stepSizeScaling reset should be designed for this purpose)
                    traits.stepSizeScaling /= 2; //TODO check if there is the necessity to multiply or reset
                    traits.nonDecreasingIterations = 0u;
                    _VV(options, "Previous iterations had no improvement. Set MY parameter to " << traits.stepSizeScaling);
                }

                if (prev_stepSize != traits.stepSize)
                    _VV(options, "The new step size for alignment " << i << " is " << traits.stepSize);

                // Update Lambda Steps using gamma of the new lines
                updateLambdaStep(traits, listUnclosedLoopMask);
            }
            saveBestAligns(traits, alignments[i], resultsSimd[i], iter);
        }

        if (foundAnOptimalAlignment)
        {
            for (std::size_t i = isOptimal.size(); i > 0; --i)
            {
                if (isOptimal[i - 1])
                {
                    goldRnaAligns.push_back(aliTraitsVect[i - 1]);
                    aliTraitsVect.erase(aliTraitsVect.begin() + i - 1);
                    erase(alignments, i - 1);
                    erase(resultsSimd, i - 1);
                    isOptimal.erase(isOptimal.begin() + i - 1);
                }
            }
        }
        if (options.verbose == 1) std::cerr << "|";
    }

    if (options.verbose > 0) std::cerr << std::endl;

    if (!empty(aliTraitsVect))
    {
        _VV(options, "Plot of the aliTraitsVect structure " << std::endl);
        plotOutput(options, aliTraitsVect);
    }

    if (!empty(goldRnaAligns))
    {
        _VV(options, "Plot of the goldRnaAligns structure " << std::endl);
        plotOutput(options, goldRnaAligns);
    }

    aliTraitsVect.insert(aliTraitsVect.end(), goldRnaAligns.begin(), goldRnaAligns.end());

    if(!empty(aliTraitsVect)) // This is a multiple alignment an the T-Coffee library must be printed
    {
        createTCoffeeLib(options, filecontents, aliTraitsVect);
    }
    return 0;
}
