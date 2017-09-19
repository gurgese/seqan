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
#include <string>
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
    res = parse(options, parser, argc, argv); // Fill the options structure with the standard and the acquired arguments
    if (res != ArgumentParser::ParseResult::PARSE_OK)  // Check the arguments
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
    SEQAN_ASSERT_EQ(length(setH), length(setV));
    SEQAN_ASSERT_EQ(length(setH), length(rnaAligns));
    _VV(options, "Number of pairwise alignments to be computed: " << length(rnaAligns));

    // A StringSet of seqan alignments is created
    StringSet<TAlign> alignsSimd;
    createSeqanAlignments(alignsSimd, setH, setV);

    // first non-structural alignment is computed
    String<TScoreValue> resultsSimd;
    _VV(options, "Start first alignment...")
    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);
    SEQAN_ASSERT_EQ(length(alignsSimd), length(rnaAligns));
    SEQAN_ASSERT_EQ(length(alignsSimd), length(resultsSimd));
    _VV(options, "\ninitial non-structural alignment (score " << resultsSimd[0] << "):\n" << alignsSimd[0]);

    // bitvector that expresses whether an alignment is finished (i.e. bound difference is very small)
    std::vector<bool> eraseV;
    eraseV.resize(length(rnaAligns), false);
    bool checkEraseV = false;

    for (TRnaAlign & ali : rnaAligns)
    {
        // apply scaling of the score matrix, according to run time parameter ssc
        //TODO check if scaling is not performed again anywhere else
        ali.structScore.score_matrix = options.laraScoreMatrix;
        ali.structScore.score_matrix.data_gap_extend /= options.sequenceScale;
        ali.structScore.score_matrix.data_gap_open   /= options.sequenceScale;
        for (unsigned j = 0; j < length(options.laraScoreMatrix.data_tab[j]); ++j)
            ali.structScore.score_matrix.data_tab[j] /= options.sequenceScale;

        // initialize memory for the fields of the alignment property data structure
        std::size_t const max_seq_size = std::max(numVertices(ali.bppGraphH.inter), numVertices(ali.bppGraphV.inter));
        std::size_t const min_seq_size = std::min(numVertices(ali.bppGraphH.inter), numVertices(ali.bppGraphV.inter));
        resize(ali.lamb, max_seq_size);           // length of longer sequence
        reserve(ali.mask, min_seq_size);           // length of shorter sequence
        resize(ali.upperBoundVect, numVertices(ali.bppGraphV.inter));
        ali.my = options.my;

        // set pointer to lambda vector
        ali.structScore.lamb = & ali.lamb;
    }

    // START FIRST ITERATION
    //#pragma omp parallel for num_threads(options.threads)
    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        // create mask of current alignment to be used for the upper/lower bound computation and the lambda update
        createMask(rnaAligns[i], alignsSimd[i]);
        std::cerr << "Created mask:";
        for (auto & mask_pair : rnaAligns[i].mask)
            std::cerr << " (" << mask_pair.first << "," << mask_pair.second << ")";
        std::cerr << std::endl;

        if (options.lowerBoundMethod == LBLEMONMWM) // The MWM is computed to fill the LowerBound
        {
            // data structure that will be passed to the lemon::MWM function to compute the full lowerBound
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(numVertices(rnaAligns[i].bppGraphH.inter));
            computeBounds(rnaAligns[i], & lowerBound4Lemon);  // upperBoundVect receives seq indices of best pairing
            computeUpperBoundScore(rnaAligns[i]); // upperBound = sum of all probability lines
            myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
            rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
            // rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
            _VV(options, "Computed maximum weighted matching using the LEMON library.");
            _VV(options, "Initial l/u bounds = " << rnaAligns[i].lowerBound << " / " << rnaAligns[i].upperBound);
        }
        else if (options.lowerBoundMethod == LBAPPROXMWM) // Approximation of MWM is computed to fill the LowerBound
        {
            computeBounds(rnaAligns[i], NULL);
            computeLowerAndUpperBoundScore(rnaAligns[i]);
        }
        else if (options.lowerBoundMethod == LBMWMTEST) // Function used to test the aproximation of MWM is computed to fill the LowerBound
        {
            //  In this branch three different methods are available for the computation: 1) the MWM approx, 2) the lemon MWM, 3) the seqan MWM <to be implemented>
            //  The approximation is used while the other structures are computed
            //  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound

            // Compute the MWM with the Lemon library
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(length(rnaAligns[i].mask));
            computeBounds(rnaAligns[i], & lowerBound4Lemon);
            computeLowerAndUpperBoundScore(rnaAligns[i]);  // also calculate GU approximation
            myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);

            // Compute the MWM with the seqan greedy MWM algorithm
            computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);

            _VVV(options, "Upper bound              = " << rnaAligns[i].upperBound);
            _VVV(options, "Lower Bound lemon primal = " << rnaAligns[i].lowerLemonBound.mwmPrimal << " \tdual = "
                      << rnaAligns[i].lowerLemonBound.mwmDual);
            _VVV(options, "Lower bound seqan greedy = " << rnaAligns[i].lowerGreedyBound);
            _VVV(options, "Lower bound approx       = " << rnaAligns[i].lowerBound);
            _VVV(options, "num edges (slm) = " << rnaAligns[i].slm);
        }
        else if(options.lowerBoundMethod == LBLINEARTIMEMWM) // using greedy algorithm
        {
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(length(rnaAligns[i].mask));
            computeBounds(rnaAligns[i], & lowerBound4Lemon);
            computeUpperBoundScore(rnaAligns[i]);
            computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
            rnaAligns[i].lowerBound = rnaAligns[i].lowerGreedyBound;
        }

        //  Compute the step size for the Lambda update
        if (rnaAligns[i].slm > 0) // if there are unpaired edges
            rnaAligns[i].stepSize = rnaAligns[i].my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);
        else
            rnaAligns[i].stepSize = 0;
        _VV(options, "The initial step size for alignment " << i << " is " << rnaAligns[i].stepSize);

        // The alignment that gives the smallest difference between upper and lower bound is saved
        saveBestAligns(rnaAligns[i], alignsSimd[i], resultsSimd[i], -1);

        if ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon))
        {
            // alignment is finished
            eraseV[i] = true;
            checkEraseV = true;
            _VV(options, "Computation for alignment " << i << " stops in iteration 0 and the bestAlignMinBounds "
                                                      << "is returned.");
        }
        else
        {
            updateLambda(rnaAligns[i]);
        }

    }
    //printRnaStructAlign(rnaAligns[0], 0);

    // Create the alignment data structure that will host the alignments with small difference between bounds
    // Move all finished alignments to goldRnaAligns, such that they are not further processed
    TRnaAlignVect goldRnaAligns;
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

    // String<TScoringSchemeStruct> alignsSimdLamb;
    // seqan::resize(alignsSimdLamb, length(alignsSimd));
    // Add struct scoring scheme pointers to each alignment cell of the alignment vector
    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        rnaAligns[i].structScore.lamb = & rnaAligns[i].lamb;
        /* rnaAligns[i].structScore.score_matrix = options.laraScoreMatrix;
        for(unsigned j = 0; j < length(options.laraScoreMatrix.data_tab[j]); ++j)
        {
            rnaAligns[i].structScore.score_matrix.data_tab[j] = rnaAligns[i].structScore.score_matrix.data_tab[j] /
                    options.sequenceScale;
        } */
    }

    // ITERATIONS OF ALIGNMENT AND MWM
    for (unsigned iter = 0; iter < options.iterations && length(alignsSimd) > 0; ++iter)
    {
        // All structural alignment is computed
        simdAlignsGlobalLocal(resultsSimd, alignsSimd, rnaAligns, options);
        checkEraseV = false;

        //#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignsSimd); ++i)
        {
            createMask(rnaAligns[i], alignsSimd[i]);
            // The MWM is computed to fill the LowerBound
            if (options.lowerBoundMethod == LBLEMONMWM)
            {
                std::pair<double, double> old_bounds{rnaAligns[i].lowerBound, rnaAligns[i].upperBound};
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(length(rnaAligns[i].mask));
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeUpperBoundScore(rnaAligns[i]);
                myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
                rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
                // rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
                if (old_bounds.first != rnaAligns[i].lowerBound || old_bounds.second != rnaAligns[i].upperBound)
                    _VV(options, "new l/u bounds: " << rnaAligns[i].lowerBound << " / " << rnaAligns[i].upperBound);
            }
            else if (options.lowerBoundMethod == LBAPPROXMWM) // Approximation of MWM is computed to fill the LowerBound
            {
                computeBounds(rnaAligns[i], NULL);
                computeLowerAndUpperBoundScore(rnaAligns[i]);
            }
            else if (options.lowerBoundMethod == LBMWMTEST) // Function used to test the aproximation of MWM is computed to fill the LowerBound
            {
                //  In this branch three different methods are available for the computation: 1) the MWM approx, 2) the lemon MWM, 3) the seqan MWM <to be implemented>
                //  The approximation is used while the other structures are computed
                //  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound

                // Compute the MWM with the Lemon library
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(length(rnaAligns[i].mask));
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeLowerAndUpperBoundScore(rnaAligns[i]);  // also calculate GU approximation
                myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);

                // Compute the MWM with the seqan greedy MWM algorithm
                computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);

                _VVV(options, "Upper bound              = " << rnaAligns[i].upperBound);
                _VVV(options, "Lower Bound lemon primal = " << rnaAligns[i].lowerLemonBound.mwmPrimal << " \tdual = "
                                                           << rnaAligns[i].lowerLemonBound.mwmDual);
                _VVV(options, "Lower bound seqan greedy = " << rnaAligns[i].lowerGreedyBound);
                _VVV(options, "Lower bound approx       = " << rnaAligns[i].lowerBound);
                _VVV(options, "num edges (slm) = " << rnaAligns[i].slm);
            }
            else if(options.lowerBoundMethod == LBLINEARTIMEMWM) // using greedy algorithm
            {
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(length(rnaAligns[i].mask));
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeUpperBoundScore(rnaAligns[i]);

                computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
                rnaAligns[i].lowerBound = rnaAligns[i].lowerGreedyBound;
                // rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
            }

            if ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon))
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
                double const prev_stepSize = rnaAligns[i].stepSize;

                //  Compute the step size for the Lambda update
                if (rnaAligns[i].slm > 0)
                    rnaAligns[i].stepSize = rnaAligns[i].my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);
                else
                    rnaAligns[i].stepSize = 0;

                // there was nothing going on in the last couple of iterations, half rnaAligns[i].my therefore
                if (rnaAligns[i].nonDecreasingIterations >= options.nonDecreasingIterations)
                {
                    // this value in case of decreasing stepsize
                    // (an opposite mechanism or a my reset should be designed for this purpose)
                    rnaAligns[i].my /= 2; //TODO check if there is the necessity to multiply or reset
                    rnaAligns[i].nonDecreasingIterations = 0u;
                    _VV(options, "Previous iterations had no improvement. Set MY parameter to " << rnaAligns[i].my);
                }

                //  Check the number of non decreasing iterations
                if (prev_stepSize < rnaAligns[i].stepSize)
                {
                    ++rnaAligns[i].nonDecreasingIterations;
                }
                else
                {
                    //TODO evaluate if the reset of this value is the right strategy wrt. the decremental solution
                    rnaAligns[i].nonDecreasingIterations = 0u;
                }

                if (prev_stepSize != rnaAligns[i].stepSize)
                    _VV(options, "The new step size for alignment " << i << " is " << rnaAligns[i].stepSize);

                updateLambda(rnaAligns[i]);
            }
            saveBestAligns(rnaAligns[i], alignsSimd[i], resultsSimd[i], iter);
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
