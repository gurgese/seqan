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
    //TODO Smollest sequences must be placed always in setV and rnaAligns[ ].bppGraphV for using NW unbunded alignment
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
        // apply scaling of the score matrix, according to run time parameter ssc
        ali.structScore.score_matrix = options.laraScoreMatrix;
        ali.structScore.score_matrix.data_gap_extend /= options.sequenceScale;
        ali.structScore.score_matrix.data_gap_open   /= options.sequenceScale;
        for (unsigned j = 0; j < length(options.laraScoreMatrix.data_tab[j]); ++j)
            ali.structScore.score_matrix.data_tab[j] /= options.sequenceScale;

        // initialize memory for the fields of the alignment property data structure
        resize (ali.lamb, std::max(numVertices(ali.bppGraphH.inter), numVertices(ali.bppGraphV.inter)));
        // TODO implement the maximum and minimum check before to come in this loop
        //resize (ali.lamb, numVertices(ali.bppGraphH.inter));
        reserve(ali.mask, std::min(numVertices(ali.bppGraphH.inter), numVertices(ali.bppGraphV.inter)));
        // This will destroy my info concerning the indexing
        //reserve(ali.mask, numVertices(ali.bppGraphV.inter));
        resize(ali.weightLineVect, numVertices(ali.bppGraphV.inter));
        ali.my = options.my;

        // Add struct scoring scheme pointers to each alignment cell of the alignment vector
        // set pointer to lambda vector
        ali.structScore.lamb = & ali.lamb;
    }

    // Create the alignment data structure that will host the alignments with small difference between bounds
    // Move all finished alignments to goldRnaAligns, such that they are not further processed
    TRnaAlignVect goldRnaAligns;

    // Receives the alignment scores
    String<TScoreValue> resultsSimd;

    _VV(options, "Start first alignment...")
    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);
//    for (unsigned i = 0; i < length(alignsSimd); ++i)
//    {
//        TRnaAlign & ali = rnaAligns[i];
//        createMask(ali, alignsSimd[i], options, true);
//        for (unsigned j = 0; j < length(ali.mask); ++j)
//        {
//            std::cerr << ali.mask[j].first << " : " << ali.mask[j].second << std::endl;
//        }
//        computeLambda(ali, options.useOppositLineUB);
//
//    }
    //_VV(options, "\nalignment in iteration " << " (score " << resultsSimd[0] << "):\n" << alignsSimd[0]);


    // ITERATIONS OF ALIGNMENT AND MWM
    for (unsigned iter = 0; iter < options.iterations && length(alignsSimd) > 0; ++iter)
    {
        checkEraseV = false;
        if (iter != 0)
            simdAlignsGlobalLocal(resultsSimd, alignsSimd, rnaAligns, options);

        _VV(options, "\nalignment in iteration " << iter << " (score " << resultsSimd[0] << "):\n" << alignsSimd[0]);

        //#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignsSimd); ++i)
        {
            TRnaAlign & ali = rnaAligns[i];
            createMask(ali, alignsSimd[i], options, true);
            std::cerr << "Created mask:";
            for (auto & mask_pair : ali.mask)
                std::cerr << " (" << mask_pair.first << "," << mask_pair.second << ")";
            std::cerr << std::endl;

            //ali.upperBound = resultsSimd[i];

            // The MWM is computed to fill the LowerBound
            if (options.lowerBoundMethod == LBLEMONMWM)
            {
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(numVertices(ali.bppGraphH.inter)); //TODO check this
                computeLowerBound(ali, & lowerBound4Lemon);
                //return 1; // I GU check until this point
                computeBounds(ali, & lowerBound4Lemon); // weightLineVect receives seq indices of best pairing
                computeUpperBoundScore(ali); // upperBound = sum of all probability lines

                ali.currentUpperBound = resultsSimd[i] + ali.upperBound;
                ali.bestUpperBound = std::min(ali.bestUpperBound, ali.currentUpperBound);

                myLemon::computeLowerBoundScore(lowerBound4Lemon, ali);
                ali.lowerBound = ali.lowerLemonBound.mwmPrimal + ali.sequenceScore;
                ali.currentLowerBound = ali.lowerBound;
                ali.bestLowerBound = std::max(ali.bestLowerBound, ali.currentLowerBound);
                // ali.slm = ali.slm - (ali.lowerLemonBound.mwmCardinality * 2);
                _VV(options, "Computed maximum weighted matching using the LEMON library.");
                _VV(options, "l/u bounds: " << ali.currentLowerBound << " / " << ali.currentUpperBound);
                _VV(options, "best l/u bounds: " << ali.bestLowerBound << " / " << ali.bestUpperBound);
            }
            else if (options.lowerBoundMethod == LBAPPROXMWM) // Approximation of MWM is computed to fill the LowerBound
            {
                computeBounds(ali, NULL);
                computeLowerAndUpperBoundScore(ali);
            }
            else if (options.lowerBoundMethod == LBMWMTEST) // Function used to test the aproximation of MWM is computed to fill the LowerBound
            {
                //  In this branch three different methods are available for the computation: 1) the MWM approx, 2) the lemon MWM, 3) the seqan MWM <to be implemented>
                //  The approximation is used while the other structures are computed
                //  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound

                // Compute the MWM with the Lemon library
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(length(ali.mask));
                computeBounds(ali, & lowerBound4Lemon);
                computeLowerAndUpperBoundScore(ali);  // also calculate GU approximation
                myLemon::computeLowerBoundScore(lowerBound4Lemon, ali);

                // Compute the MWM with the seqan greedy MWM algorithm
                computeLowerBoundGreedy(lowerBound4Lemon, ali);

                _VVV(options, "Upper bound              = " << ali.upperBound);
                _VVV(options, "Lower Bound lemon primal = " << ali.lowerLemonBound.mwmPrimal << " \tdual = "
                                                            << ali.lowerLemonBound.mwmDual);
                _VVV(options, "Lower bound seqan greedy = " << ali.lowerGreedyBound);
                _VVV(options, "Lower bound approx       = " << ali.lowerBound);
                _VVV(options, "num edges (slm) = " << ali.slm);
            }
            else if(options.lowerBoundMethod == LBLINEARTIMEMWM) // using greedy algorithm
            {
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(length(ali.mask));
                computeBounds(ali, & lowerBound4Lemon);
                computeUpperBoundScore(ali);

                computeLowerBoundGreedy(lowerBound4Lemon, ali);
                ali.lowerBound = ali.lowerGreedyBound;
                // ali.slm = ali.slm - (ali.lowerLemonBound.mwmCardinality * 2);
            }

            // TODO move saveBestAligns call here?
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
                // save previous step size
                double const prev_stepSize = ali.stepSize;

                //  Compute the step size for the Lambda update
                if (ali.slm > 0)
                    ali.stepSize = ali.my * ((ali.bestUpperBound - ali.bestLowerBound) / ali.slm);
                else
                    ali.stepSize = 0;

                if (iter == 0)
                {
                    _VV(options, "The initial step size for alignment " << i << " is " << ali.stepSize);
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

                //  Check the number of non decreasing iterations
                if (prev_stepSize < ali.stepSize)
                {
                    ++ali.nonDecreasingIterations;
                }
                else
                {
                    //TODO evaluate if the reset of this value is the right strategy wrt. the decremental solution
                    ali.nonDecreasingIterations = 0u;
                }

                if (prev_stepSize != ali.stepSize)
                    _VV(options, "The new step size for alignment " << i << " is " << ali.stepSize);

                updateLambda(ali);
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
