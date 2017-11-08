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
    // Argument Parser
    if (argc == 1)
    {
        std::cout << "type " << argv[0] << " --help to get the parameters table (-i option is mandatory)" << std::endl;
        return 1;
    }

    ArgumentParser parser;
    Options options;
    setupArgumentParser(parser, options);
    ArgumentParser::ParseResult res;
    res = parse(options, parser, argc, argv); // Fill the options structure
    if (res != ArgumentParser::ParseResult::PARSE_OK)
        return res == ArgumentParser::ParseResult::PARSE_ERROR ? 1 : 0;

    // Read input files
    RnaStructContents filecontents1;
    _readRnaInputFile(filecontents1, options.inFile, options);
    RnaStructContents filecontents2;
    _readRnaInputFile(filecontents2, options.inFileRef, options);

    _VV(options, "Read " << length(filecontents1.records) << " and "
                         << length(filecontents2.records) << " records from input files.");

    // Add the weight interaction edges vector map in the data structure using Vienna package
    bppInteractionGraphBuild(filecontents1.records, options);
    bppInteractionGraphBuild(filecontents2.records, options);
    _VV(options, getEbpseqString(filecontents1) << getEbpseqString(filecontents2));

    // CREATE PAIRWISE ALIGNMENTS
    // If one input file is given, then build unique pairs of the input sequences.
    // If two input files are given, then build the cross product of all sequences from the first file
    //     with all sequences from the second file.
    // Of each pair the first sequence is stored in setH and the second sequence is stored in setV.
    RnaSeqSet setH;
    RnaSeqSet setV;
    TRnaAlignVect rnaAligns;
    crossproduct(setH, setV, rnaAligns, filecontents1.records, filecontents2.records);
    //TODO Smallest sequences must be placed always in setV and rnaAligns[ ].bppGraphV for using NW unbounded alignment
    SEQAN_ASSERT_EQ(length(setH), length(setV));
    SEQAN_ASSERT_EQ(length(setH), length(rnaAligns));
    _VV(options, "Number of pairwise alignments to be computed: " << length(rnaAligns));

    // A StringSet of seqan alignments is created
    StringSet<TAlign> alignsSimd;
    createSeqanAlignments(alignsSimd, setH, setV);
    SEQAN_ASSERT_EQ(length(alignsSimd), length(rnaAligns));

    // Bitvector that expresses whether an alignment is finished (i.e. bound difference is very small)
    std::vector<bool> eraseV;
    eraseV.resize(length(rnaAligns), false);
    bool checkEraseV = false;

    for (TRnaAlign & ali : rnaAligns)
    {
        // the lambda structs are referenced by the index of the first sequence
        resize(ali.lamb, numVertices(ali.bppGraphH.inter));
        // the weight lines are referenced by the index of the second sequence
        resize(ali.weightLineVect, numVertices(ali.bppGraphV.inter));
        // allocate capacity for mask: the length of the shorter sequence is the maximum number of expressed lines
        reserve(ali.mask, std::min(numVertices(ali.bppGraphH.inter), numVertices(ali.bppGraphV.inter)));

        ali.structScore.score_matrix = options.laraScoreMatrix;
        ali.structScore.lamb = & ali.lamb;
        ali.my = options.my;
    }

    // Create the alignment data structure that will host the alignments with small difference between bounds
    // Move all finished alignments to goldRnaAligns, such that they are not further processed
    TRnaAlignVect goldRnaAligns;

    // Receives the alignment scores
    String<TScoreValue> resultsSimd;

    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);
    _VV(options, "\ninitial sequence alignment (score " << resultsSimd[0] << "):\n" << alignsSimd[0]);

    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        TRnaAlign & ali = rnaAligns[i];
        createMask(ali, alignsSimd[i], options);
//        for (unsigned j = 0; j < length(ali.mask); ++j)
//        {
//            std::cerr << ali.mask[j].first << " : " << ali.mask[j].second << std::endl;
//        }

        createNewLambdaLines(ali, options.useOppositLineUB);
    }


    // ITERATIONS OF ALIGNMENT AND MWM
    for (unsigned iter = 0; iter < options.iterations && length(alignsSimd) > 0; ++iter)
    {
        checkEraseV = false;
        //if (iter != 0)
        simdAlignsGlobalLocal(resultsSimd, alignsSimd, rnaAligns, options);
        _VV(options, "\nalignment in iteration " << iter << " (score " << resultsSimd[0] << "):\n" << alignsSimd[0]);

        //#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignsSimd); ++i)
        {
            TRnaAlign & ali = rnaAligns[i];
            bool changedMask = createMask(ali, alignsSimd[i], options);

            ali.upperBound = resultsSimd[i];
            ali.bestUpperBound = std::min(ali.bestUpperBound, ali.upperBound);
//            std::cerr << "Saved Upper Bound (alignment Score)" << std::endl;

            if (changedMask || i == 0)
            {

                createNewLambdaLines(ali, options.useOppositLineUB);
//                std::cerr << "Included in Lambda vector the new alignment lines" << std::endl;

                // The MWM is computed to fill the LowerBound
                if (options.lowerBoundMethod == LBLEMONMWM) {
                    TMapVect lowerBound4Lemon;
                    lowerBound4Lemon.resize(numVertices(ali.bppGraphH.inter)); //TODO check this
                    //std::cout << alignsSimd[i] << std::endl;
                    computeLowerBound(ali, & lowerBound4Lemon);
//                computeBounds(ali, & lowerBound4Lemon); // weightLineVect receives seq indices of best pairing
//                computeUpperBoundScore(ali); // upperBound = sum of all probability lines
                    myLemon::computeLowerBoundScore(lowerBound4Lemon, ali);
                    ali.lowerBound = ali.lowerLemonBound.mwmPrimal + ali.sequenceScore;
                    // ali.slm = ali.slm - (ali.lowerLemonBound.mwmCardinality * 2);
                    _VV(options, "Computed maximum weighted matching using the LEMON library.");
                    ali.bestLowerBound = std::max(ali.bestLowerBound, ali.lowerBound);
                }
                else if (options.lowerBoundMethod == LBAPPROXMWM)
                {   // TODO Verify the procedure! Approximation of MWM is computed to fill the LowerBound
                    computeBounds(ali, NULL);  // TODO verify this call and reimplement the approximation if better than gready
                } else if (options.lowerBoundMethod == LBLINEARTIMEMWM) // TODO Verify the procedure! using greedy algorithm
                {
                    TMapVect lowerBound4Lemon;
                    lowerBound4Lemon.resize(length(ali.mask));
                    computeBounds(ali, &lowerBound4Lemon);
                    computeLowerBoundGreedy(lowerBound4Lemon, ali);
                    ali.lowerBound = ali.lowerGreedyBound;
                    // ali.slm = ali.slm - (ali.lowerLemonBound.mwmCardinality * 2);
                }
            }
            _VV(options,"best l/u bounds: " << ali.bestLowerBound << " / " << ali.bestUpperBound);

            if ((ali.bestUpperBound - ali.bestLowerBound < options.epsilon))
            {
                // alignment is finished
                eraseV[i] = true;
                checkEraseV = true;
                _VV(options, "Computation for alignment " << i << " stops in iteration "
                                                          << iter << " and the bestAlignMinBounds is returned.");
            }
            else
            {
                seqan::String<std::pair <unsigned, unsigned> > listUnclosedLoopMask;
                // Compute slm factor
                computeSlm(ali, listUnclosedLoopMask);

                // save previous step size
                double const prev_stepSize = ali.stepSize;
                //  Compute the step size for the Lambda update
                if (ali.slm > 0)
                    ali.stepSize = ali.my * ((ali.bestUpperBound - ali.bestLowerBound) / ali.slm);
                else
                    ali.stepSize = 0;


                //  Check the number of non decreasing iterations
                if (prev_stepSize <= ali.stepSize)
                {
                    ++ali.nonDecreasingIterations;
                    _VV(options, "nonDecreasingIterations = " << ali.nonDecreasingIterations);
                }
                else
                {
                    //TODO evaluate if the reset of this value is the right strategy wrt. the decremental solution
                    ali.nonDecreasingIterations = 0u;
                }

                // there was nothing going on in the last couple of iterations, half ali.my therefore
                if (ali.nonDecreasingIterations >= options.nonDecreasingIterations)
                {
                    // this value in case of decreasing stepsize
                    // (an opposite mechanism or a my reset should be designed for this purpose)
                    ali.my /= 2; //TODO check if there is the necessity to multiply or reset
                    ali.nonDecreasingIterations = 0u;
                    _VV(options, "Previous iterations had no improvement. Set MY parameter to " << ali.my);
                }

                if (prev_stepSize != ali.stepSize)
                    _VV(options, "The new step size for alignment " << i << " is " << ali.stepSize);

                // Update Lambda Steps using gamma of the new lines
                updateLambdaStep(ali, listUnclosedLoopMask);
            }
            saveBestAligns(ali, alignsSimd[i], resultsSimd[i], iter);
        }

        if (checkEraseV)
        {
            for (std::size_t i = eraseV.size(); i > 0; --i)
            {
                if (eraseV[i - 1])
                {
                    goldRnaAligns.push_back(rnaAligns[i - 1]);
                    rnaAligns.erase(rnaAligns.begin() + i - 1);
                    erase(alignsSimd, i - 1);
                    erase(resultsSimd, i - 1);
                    eraseV.erase(eraseV.begin() + i - 1);
                }
            }
        }
        if (options.verbose == 1) std::cerr << "|";
    }

    if (options.verbose > 0) std::cerr << std::endl;

    if (!empty(rnaAligns))
    {
        _VV(options, "Plot of the rnaAligns structure " << std::endl);
        plotOutput(options, rnaAligns);
    }

    if (!empty(goldRnaAligns))
    {
        _VV(options, "Plot of the goldRnaAligns structure " << std::endl);
        plotOutput(options, goldRnaAligns);
    }

    rnaAligns.insert(rnaAligns.end(), goldRnaAligns.begin(), goldRnaAligns.end());

    if(!empty(rnaAligns)) // This is a multiple alignment an the T-Coffee library must be printed
    {
        createTCoffeeLib(options, filecontents1, filecontents2, rnaAligns);
    }
    return 0;
}
