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
#include <omp.h>
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
    setScoreMatrix(options);

    // Read input files.
    // If one input file is given, then build unique pairs of the input sequences.
    // If two input files are given, then build the cross product of all sequences from the first file
    //     with all sequences from the second file.
    RnaStructContentsPair filecontents{};
    readRnaFile(filecontents.first, options.inFile, options);
    readRnaFile(filecontents.second, options.inFileRef, options);
    _VV(options, "Read " << length(filecontents.first.records) << " and "
                         << length(filecontents.second.records) << " records from input files.");

    if (length(options.dotplotFile) == length(filecontents.first.records) + length(filecontents.second.records))
    {
        _V(options, "Using given dotplot files to extract the base pair probabilities.");
        readDotplotFiles(filecontents, options.dotplotFile);
    }
    else
    {
        // If not present, compute the weighted interaction edges using ViennaRNA functions.
        computeMissingInteractions(filecontents.first.records, options);
        computeMissingInteractions(filecontents.second.records, options);
    }
    _VVV(options, getEbpseqString(filecontents.first) << getEbpseqString(filecontents.second));


    // A StringSet of SeqAn-style alignments and a vector of alignment traits are created.
    StringSet<RnaAlignment> alignments;
    RnaAlignmentTraitsVector alignmentTraits;
    prepareAlignmentVector(alignments, alignmentTraits, filecontents);
    _VV(options, "Number of pairwise alignments to be computed: " << length(alignments));

    // Initialise the alignment traits.
    for (RnaAlignmentTraits & traits : alignmentTraits)
    {
        // The lambda structs are referenced by the index of the first sequence.
        resize(traits.interactions, numVertices(traits.bppGraphH.inter));
        // Allocate capacity for lines: the length of the shorter sequence is the maximum number of expressed lines.
        reserve(traits.lines, std::min(numVertices(traits.bppGraphH.inter), numVertices(traits.bppGraphV.inter)));

        traits.structureScore.matrix = options.laraScoreMatrix;
        traits.structureScore.interactions = & traits.interactions;
        traits.stepSizeScaling = options.stepSizeScaling;
    }

    // Create the alignment data structure that will host the alignments with small difference between bounds.
    // Move all finished alignments to optimalAlignTraits, such that they are not further processed.
    RnaAlignmentTraitsVector optimalAlignTraits;

    // Bitvector that expresses whether an alignment is finished (i.e. bound difference is very small)
    std::vector<bool> isOptimal;
    isOptimal.resize(length(alignments), false);

    // Calculate the initial alignment scores.
    String<double> upperBoundScores;
    initialAlignment(upperBoundScores, alignments, options);
    _VV(options, "\ninitial sequence alignment (score " << upperBoundScores[0] << "):\n" << alignments[0]);

    // Evaluate the lines and interactions.
    for (unsigned idx = 0u; idx < length(alignments); ++idx)
    {
        evaluateLines(alignmentTraits[idx], alignments[idx], options);
        evaluateInteractions(alignmentTraits[idx], 0u);
    }

    // ITERATIONS OF ALIGNMENT AND MWM
    for (unsigned iter = 0u; iter < options.iterations && length(alignments) > 0u; ++iter)
    {
        bool foundAnOptimalAlignment = false;
        structuralAlignment(upperBoundScores, alignments, alignmentTraits, options);

        #pragma omp parallel for num_threads(options.threads)
        for (unsigned idx = 0; idx < length(alignments); ++idx)
        {
            RnaAlignmentTraits & traits = alignmentTraits[idx];
            bool const changedLines = evaluateLines(traits, alignments[idx], options);

            // Set upper bound (i.e. the relaxed solution = Lagrangian dual).
            traits.upperBound = upperBoundScores[idx];
//            traits.bestUpperBound = std::min(traits.bestUpperBound, traits.upperBound); // this value must be changed in agreement with the lower bound

            // If lines have changed the lower bound must be calculated.
            if (changedLines || iter == 0)
            {
                _VV(options, "\nalignment " << idx << " in iteration " << iter << " (score " << upperBoundScores[0]
                                            << "):\n" << alignments[idx]);
                evaluateInteractions(traits, iter);

                InteractionScoreMap validInteractionScores;
                // The scores are referenced by the index of the first sequence.
                validInteractionScores.resize(numVertices(traits.bppGraphH.inter));
                prepareLowerBoundScores(validInteractionScores, traits);

                switch (options.lowerBoundMethod)
                {
                    case MWM_LEMON:  traits.lowerBound = myLemon::computeLowerBoundScore(validInteractionScores);
                                     break;
                    case MWM_GREEDY: traits.lowerBound = computeLowerBoundGreedy(validInteractionScores);
                                     break;
                    case MWM_SIMPLE: // computeBounds(traits, NULL);
                    default:         std::cerr << "Lower Bound method not implemented." << std::endl;
                                     exit(1);
                }
                _VVV(options,iter << "\tLow bound component MWM/SeqAlignScore for sequences (" << traits.sequenceIndices.first
                                                              << ","  << traits.sequenceIndices.second << "): "
                                                              << traits.lowerBound << " / " << traits.sequenceScore);
                traits.lowerBound += traits.sequenceScore;
            }
            // Save best alignment scores for smallest Epslon
            if( (traits.upperBound - traits.lowerBound) < (traits.bestUpperBound - traits.bestLowerBound) )
            {
                traits.bestUpperBound = traits.upperBound;
                traits.bestLowerBound = traits.lowerBound;
            }
            // Save best alignment scores for maximum lowerBound
            if( traits.lowerBound > traits.bestLowerBoundMaxLow)
            {
                traits.bestUpperBoundMaxLow = traits.upperBound;
                traits.bestLowerBoundMaxLow = traits.lowerBound;
            }
            // Save best alignment scores for minimum upperBound
            if( traits.upperBound < traits.bestUpperBoundMinUp) // TODO verify if there is the possibility to enter here without updating the lowerBound
            {
                traits.bestUpperBoundMinUp = traits.upperBound;
                traits.bestLowerBoundMinUp = traits.lowerBound;
            }
//            traits.bestUpperBound = std::min(traits.bestUpperBound, traits.upperBound);
//            traits.bestLowerBound = std::max(traits.bestLowerBound, traits.lowerBound);
            _VVV(options,iter << "\tCurrent l/u bounds for sequences (" << traits.sequenceIndices.first << ","
                                                          << traits.sequenceIndices.second << "): "
                                                          << traits.lowerBound << " / " << traits.upperBound);
            _VVV(options,iter << "\tBest l/u bounds with maximum LowerBound for sequences (" << traits.sequenceIndices.first << ","
                              << traits.sequenceIndices.second << "): "
                              << traits.bestLowerBoundMaxLow << " / " << traits.bestUpperBoundMaxLow);
            _VVV(options,iter << "\tBest l/u bounds with minimum UpperBound for sequences (" << traits.sequenceIndices.first << ","
                             << traits.sequenceIndices.second << "): "
                             << traits.bestLowerBoundMinUp << " / " << traits.bestUpperBoundMinUp);
            _VV(options,iter << "\tBest l/u bounds with smallest epslon for sequences (" << traits.sequenceIndices.first << ","
                             << traits.sequenceIndices.second << "): "
                             << traits.bestLowerBound << " / " << traits.bestUpperBound);

            SEQAN_ASSERT_LEQ(traits.bestLowerBound, traits.bestUpperBound);

            if (traits.bestUpperBound - traits.bestLowerBound < options.epsilon)
            {
                // This alignment is already optimal.
                isOptimal[idx] = true;
                foundAnOptimalAlignment = true;
                _VV(options, "Computation for alignment " << idx << " stops in iteration " << iter << ".");
            }
            else
            {
                // Compute the number of active subgradients (s_lm).
                seqan::String<PositionPair> unclosedLoops;
                computeNumberOfSubgradients(traits, unclosedLoops);

                //  Compute the step size for the subgradient update.
                double const previousStepSize = traits.stepSize;
                if (traits.numberOfSubgradients > 0)
                    traits.stepSize = traits.stepSizeScaling * ((traits.bestUpperBound - traits.bestLowerBound)
                                                                / traits.numberOfSubgradients);
                else
                    traits.stepSize = 0;

                //  Check the number of non decreasing iterations
                if (previousStepSize <= traits.stepSize)
                {
                    ++traits.nonDecreasingIterations;
                    _VVV(options, "nonDecreasingIterations = " << traits.nonDecreasingIterations);
                }
                else
                {
                    traits.nonDecreasingIterations = 0u;
                }

                // There was no improvement in the last couple of iterations, therefore half traits.stepSizeScaling.
                if (traits.nonDecreasingIterations >= options.nonDecreasingIterations)
                {
                    // this value in case of decreasing stepsize
                    // (an opposite mechanism or a stepSizeScaling reset should be designed for this purpose)
                    traits.stepSizeScaling /= 2;
                    traits.nonDecreasingIterations = 0u;
                    _VV(options, "Previous iterations had no improvement. " << "Set step size scaling parameter to "
                                                                            << traits.stepSizeScaling);
                }

                // Update subgradients using the new unclosed interactions.
                updateLambdaValues(traits, unclosedLoops);
            }
            saveBestAligns(traits, alignments[idx], upperBoundScores[idx], iter);
        }

        if (foundAnOptimalAlignment)
        {
            for (std::size_t idx = isOptimal.size(); idx > 0; --idx)
            {
                if (isOptimal[idx - 1])
                {
                    optimalAlignTraits.push_back(alignmentTraits[idx - 1]);
                    alignmentTraits.erase(alignmentTraits.begin() + idx - 1);
                    erase(alignments, idx - 1);
                    erase(upperBoundScores, idx - 1);
                    isOptimal.erase(isOptimal.begin() + idx - 1);
                }
            }
            // renew pointers
            for (RnaAlignmentTraits & traits : alignmentTraits)
                traits.structureScore.interactions = & traits.interactions;
        }
        // Progress bar.
        if (options.verbose == 1) std::cerr << "|";
    }

    if (options.verbose > 0) std::cerr << std::endl;

    if (!empty(alignmentTraits))
    {
        _VV(options, "Plot of the alignmentTraits structure " << std::endl);
        plotAlignments(options, alignmentTraits);
    }

    if (!empty(optimalAlignTraits))
    {
        _VV(options, "Plot of the optimalAlignTraits structure " << std::endl);
        plotAlignments(options, optimalAlignTraits);
    }

    alignmentTraits.insert(alignmentTraits.end(), optimalAlignTraits.begin(), optimalAlignTraits.end());

    // This is a multiple alignment and the T-Coffee library must be printed
    if (endsWith(options.outFile, "fa") || endsWith(options.outFile, "fasta"))
        createFastaAlignmentFile(options, alignmentTraits);
    else
        createTCoffeeLib(options, filecontents, alignmentTraits);

    return 0;
}
