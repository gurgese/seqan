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

    _V(options, "Read " << length(filecontents1.records) << " and "
                        << length(filecontents2.records) << " records from input files.");

    // Add the weight interaction edges vector map in the data structure using Vienna package
    bppInteractionGraphBuild(filecontents1.records, options);
    bppInteractionGraphBuild(filecontents2.records, options);
    std::cout << filecontents1.records.front().bppMatrGraphs[0].inter << std::endl;

    // create pairwise alignments
    RnaSeqSet setH;
    RnaSeqSet setV;
    TRnaAlignVect rnaAligns;
    crossproduct(setH, setV, rnaAligns, filecontents1.records, filecontents2.records);

//  Create the alignment data structure that will host the alignments with small difference between upper and lower bound
    TRnaAlignVect goldRnaAligns;
    std::vector<bool> eraseV;
    bool checkEraseV = false;

    StringSet<TAlign> alignsSimd;
    String<TScoreValue> resultsSimd;
    // simd vector is created
    createSimdAligns(alignsSimd, setH, setV);
    eraseV.assign(length(alignsSimd), false);

// timer start
    std::clock_t begin = std::clock();
    std::chrono::steady_clock::time_point beginChrono = std::chrono::steady_clock::now();

// first non-structural alignment is computed
    firstSimdAlignsGlobalLocal(resultsSimd, alignsSimd, options);

    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        rnaAligns[i].structScore.score_matrix = options.laraScoreMatrix;
        for(unsigned j = 0; j < length(options.laraScoreMatrix.data_tab[j]); ++j)
        {
            rnaAligns[i].structScore.score_matrix.data_tab[j] = rnaAligns[i].structScore.score_matrix.data_tab[j] /
                                                                options.sequenceScale;
            rnaAligns[i].structScore.score_matrix.data_gap_open = rnaAligns[i].structScore.score_matrix.data_gap_open /
                                                                options.sequenceScale;
            rnaAligns[i].structScore.score_matrix.data_gap_extend = rnaAligns[i].structScore.score_matrix.data_gap_extend /
                                                                    options.sequenceScale;
//TODO sequenceScale can be substituted from a runtime computed parameter that consider the identity of the sequences or other aspects
        }
    }

    double mwmtime = 0.0;
    double lemtime = 0.0;
    double boutime = 0.0;

//#pragma omp parallel for num_threads(options.threads)
    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        resize(rnaAligns[i].lamb, length(setH[i]));  // length of longer sequence
        resize(rnaAligns[i].mask, length(setV[i]));  // length of shorter sequence
        resize(rnaAligns[i].upperBoundVect, length(setV[i]));  // length of shorter sequence
        rnaAligns[i].my = options.my;

        SEQAN_ASSERT_LEQ(length(setV[i]), length(setH[i]));

// Save the best alignments that give the absolute maximum score
//        saveBestAlign(rnaAligns[i], alignsSimd[i], resultsSimd[i]);
// Create the mask of the current alignment to be used for the upper, lower bound computation and the lamb update
        maskCreator(rnaAligns[i], alignsSimd[i]);

        if(options.lowerBoundMethod == LBLEMONMWM) // The MWM is computed to fill the LowerBound
        {
            _VV(options, "using Lemon MWM");
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
            computeBounds(rnaAligns[i], & lowerBound4Lemon);  // upperBoundVect receives seq indices of best pairing
            computeUpperBoundScore(rnaAligns[i]);
// Compute the MWM with the Lemon library
            myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
            rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
//            rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
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
            lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
            std::clock_t clstart = std::clock();
            computeBounds(rnaAligns[i], & lowerBound4Lemon);
            computeLowerAndUpperBoundScore(rnaAligns[i]);  // also calculate GU approximation
            boutime += double(std::clock() - clstart) / CLOCKS_PER_SEC;
            clstart = std::clock();
            myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
            lemtime += double(std::clock() - clstart) / CLOCKS_PER_SEC;

            // Compute the MWM with the seqan greedy MWM algorithm
            clstart = std::clock();
            computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
            mwmtime += double(std::clock() - clstart) / CLOCKS_PER_SEC;

            _VV(options, "Upper bound              = " << rnaAligns[i].upperBound);

            _VV(options, "Lower Bound lemon primal = " << rnaAligns[i].lowerLemonBound.mwmPrimal << " \tdual = "
                      << rnaAligns[i].lowerLemonBound.mwmDual);
            _VV(options, "Lower bound seqan greedy = " << rnaAligns[i].lowerGreedyBound);
            _VV(options, "Lower bound approx       = " << rnaAligns[i].lowerBound);
            _VV(options, "num edges (slm) = " << rnaAligns[i].slm);
        }
        else if(options.lowerBoundMethod == LBLINEARTIMEMWM) // using greedy algorithm
        {
            TMapVect lowerBound4Lemon;
            lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
            computeBounds(rnaAligns[i], & lowerBound4Lemon);
            computeUpperBoundScore(rnaAligns[i]);
// Compute the MWM
            computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
            rnaAligns[i].lowerBound = rnaAligns[i].lowerGreedyBound;
        }

        //  Compute the step size for the Lambda update
        if(rnaAligns[i].slm > 0)
            rnaAligns[i].stepSize = rnaAligns[i].my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);
        else
            rnaAligns[i].stepSize = 0;

        int index = -1;
// The alignemnt that give the smallest difference between up and low bound should be saved
        saveBestAligns(rnaAligns[i], alignsSimd[i], resultsSimd[i], index);
//        saveBestAlignMinBound(rnaAligns[i], alignsSimd[i], resultsSimd[i], index);
//        saveBestAlignScore(rnaAligns[i], alignsSimd[i], resultsSimd[i], index);

        if ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon))
        {
            _VV(options, "Computation for this alignment should stopped and the bestAlignMinBounds should be returned "
                      "upper bound = " << rnaAligns[i].upperBound << " lower bound = " << rnaAligns[i].lowerBound);

/*
            std::cout << rnaAligns[i].idBppSeqH << "/" << rnaAligns[i].idBppSeqV <<
            " upper bound = " << rnaAligns[i].upperBound << " lower bound = " << rnaAligns[i].lowerBound << std::endl;
            if (rnaAligns[i].idBppSeqH == 3 && rnaAligns[i].idBppSeqV == 4)
            {
                std::cout << rnaAligns[i].forMinBound.bestAlign << std::endl;
                std::cout << "Graph H \n" << rnaAligns[i].bppGraphH.inter << std::endl;
                std::cout << "Graph V \n" << rnaAligns[i].bppGraphV.inter << std::endl;
                std::cout << "3 - 312 / 4 - 310 " << rnaAligns[i].lamb[312].map[310].step << " / "
                          << rnaAligns[i].lamb[312].map[310].maxProbScoreLine
                          << " : " << length(rnaAligns[i].lamb[310]) << std::endl;
                std::cout << "length 308 " << length(rnaAligns[i].lamb[308]) << " / "
                          << std::endl;

                std::cout << "length fixed " << length(filecontents1.records[rnaAligns[i].idBppSeqH].fixedGraphs[0]) << std::endl;

                std::cout << outDegree(filecontents1.records[3].bppMatrGraphs[0].inter, 308) << std::endl;

                String<unsigned int> vect;
                getVertexAdjacencyVector(vect, filecontents1.records[3].bppMatrGraphs[0].inter, 300);

                std::cout << seqan::length(vect) << std::endl;
                std::cout << degree(filecontents1.records[3].bppMatrGraphs[0].inter, 312) << std::endl;
                std::cout << degree(filecontents1.records[4].bppMatrGraphs[0].inter, 310) << std::endl;

//                TEdgeStump * currentOut1 = filecontents1.records[3].bppMatrGraphs[0].inter.data_vertex[312];
//                TEdgeStump * currentOut2 = filecontents1.records[3].bppMatrGraphs[0].inter.data_vertex[309];
//                std::cout << "312 = " << currentOut1 << " 309 = " << currentOut2 << std::endl;
            }
*/

//            eraseVect.push_back(i);
            eraseV[i] = true;
            checkEraseV = true;
        }
        else
        {

            _VVV(options, "\nThe step size to be used for Lambda for alignment " << i << " in iteration 0 is " << rnaAligns[i].stepSize);

            updateLambda(rnaAligns[i]);
        }

    }
    printRnaStructAlign(rnaAligns[0], 0);

    if (checkEraseV)
    {
        for (int i = eraseV.size() - 1; i >= 0; --i)
        {
            if(eraseV[i])
            {
                goldRnaAligns.push_back(rnaAligns[i]);
                rnaAligns.erase(rnaAligns.begin() + i);
                erase(alignsSimd, i);
                erase(resultsSimd, i);
                eraseV.erase(eraseV.begin() + i);
            }
        }
    }

//    String<TScoringSchemeStruct> alignsSimdLamb;
//    seqan::resize(alignsSimdLamb, length(alignsSimd));
// Add struct scoring scheme pointers to each alignment cell of the alignment vector
    for (unsigned i = 0; i < length(alignsSimd); ++i)
    {
        rnaAligns[i].structScore.lamb = & rnaAligns[i].lamb;
        rnaAligns[i].structScore.score_matrix = options.laraScoreMatrix;
        for(unsigned j = 0; j < length(options.laraScoreMatrix.data_tab[j]); ++j)
        {
            rnaAligns[i].structScore.score_matrix.data_tab[j] = rnaAligns[i].structScore.score_matrix.data_tab[j] /
                    options.sequenceScale;
//TODO sequenceScale can be substituted from a runtime computed parameter that consider the identity of the sequences or other aspects
        }
    }

    for (unsigned x = 0; x < options.iterations && length(alignsSimd) > 0; ++x)
    {
        // All structural alignment is computed
        simdAlignsGlobalLocal(resultsSimd, alignsSimd, rnaAligns, options);
        checkEraseV = false;
//#pragma omp parallel for num_threads(options.threads)
        for (unsigned i = 0; i < length(alignsSimd); ++i) // TODO replace this function with the SIMD implementation for execute in PARALLEL
        {

            if(options.lowerBoundMethod == LBLEMONMWM) // The MWM is computed to fill the LowerBound
            {
//  Define the datastructure that will be passed to the lemon::MWM function to compute the full lowerBound
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeUpperBoundScore(rnaAligns[i]);
// Compute the MWM with the Lemon library
                myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
                rnaAligns[i].lowerBound = rnaAligns[i].lowerLemonBound.mwmPrimal;
//                rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
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
                lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
                std::clock_t clstart = std::clock();
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeLowerAndUpperBoundScore(rnaAligns[i]);  // also calculate GU approximation
                boutime += double(std::clock() - clstart) / CLOCKS_PER_SEC;
                clstart = std::clock();
                myLemon::computeLowerBoundScore(lowerBound4Lemon, rnaAligns[i]);
                lemtime += double(std::clock() - clstart) / CLOCKS_PER_SEC;

                // Compute the MWM with the seqan greedy MWM algorithm
                clstart = std::clock();
                computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
                mwmtime += double(std::clock() - clstart) / CLOCKS_PER_SEC;

                _VV(options, "Upper bound              = " << rnaAligns[i].upperBound);

                _VV(options, "Lower Bound lemon primal = " << rnaAligns[i].lowerLemonBound.mwmPrimal << " \tdual = "
                                                           << rnaAligns[i].lowerLemonBound.mwmDual);
                _VV(options, "Lower bound seqan greedy = " << rnaAligns[i].lowerGreedyBound);
                _VV(options, "Lower bound approx       = " << rnaAligns[i].lowerBound);
                _VV(options, "num edges (slm) = " << rnaAligns[i].slm);
            }
            else if(options.lowerBoundMethod == LBLINEARTIMEMWM) // using greedy algorithm
            {
                TMapVect lowerBound4Lemon;
                lowerBound4Lemon.resize(rnaAligns[i].maskIndex);
                computeBounds(rnaAligns[i], & lowerBound4Lemon);
                computeUpperBoundScore(rnaAligns[i]);

                computeLowerBoundGreedy(lowerBound4Lemon, rnaAligns[i]);
                rnaAligns[i].lowerBound = rnaAligns[i].lowerGreedyBound;
//                rnaAligns[i].slm = rnaAligns[i].slm - (rnaAligns[i].lowerLemonBound.mwmCardinality * 2);
            }

            //  Compute the step size for the Lambda update
            double stepSize;
            if(rnaAligns[i].slm > 0)
                stepSize = rnaAligns[i].my * ((rnaAligns[i].upperBound - rnaAligns[i].lowerBound) / rnaAligns[i].slm);
            else
                stepSize = 0;

            if (rnaAligns[i].upperBound - rnaAligns[i].lowerBound < options.epsilon)
            {
                _VV(options, "Computation for this alignment should stopped and the bestAlignMinBounds should be returned "
                    "upper bound = " << rnaAligns[i].upperBound << " lower bound = " << rnaAligns[i].lowerBound);
//                eraseVect.push_back(i);
                eraseV[i] = true;
                checkEraseV = true;
            }
            else
            {
                // there was nothing going on in the last couple of iterations, half rnaAligns[i].my therefore
                if (rnaAligns[i].nonDecreasingIterations == options.nonDecreasingIterations)
                {
                    rnaAligns[i].my = rnaAligns[i].my/2; //TODO check if there is the necessity to multiply or reset
                    // this value in case of decreasing stepsize (an opposite mechanism or a my reset should be designed for this purpose)
                }

                //  Check the number of non decreasing iterations
                if(rnaAligns[i].stepSize < stepSize)
                {
                    ++rnaAligns[i].nonDecreasingIterations;
                }
                else
                {
                    rnaAligns[i].nonDecreasingIterations = 0u; //TODO evaluate if the reset of this value is the right strategy with respect to the decremental solution
                }

                // Assign the new stepSize to for the Lambda update
                rnaAligns[i].stepSize = stepSize;

                _VVV(options, "\nThe step size to be used for Lambda for alignment " << i << " in iteration" << x << " is " << rnaAligns[i].stepSize);

                updateLambda(rnaAligns[i]);
            }
            // The alignemnt that give the smallest difference between up and low bound should be saved
            saveBestAligns(rnaAligns[i], alignsSimd[i], resultsSimd[i], x);
//            saveBestAlignMinBound(rnaAligns[i], alignsSimd[i], resultsSimd[i], x);
        }
        if (x == 12 || x == 13 || x == 14 ||x == 99 || x == 499)
        {
            printRnaStructAlign(rnaAligns[0], x);

            //--------------------
            std::cout << x << "\talign row1: " << length(row(rnaAligns[0].forMinBound.bestAlign, 0)) << " ";
            for (auto c : row(rnaAligns[0].forMinBound.bestAlign, 0))
                std::cout << c;
            std::cout << "\tscore = " << resultsSimd[0] << std::endl << "\talign row2: " << length(row(rnaAligns[0].forMinBound.bestAlign, 1)) << " ";
            for (auto c : row(rnaAligns[0].forMinBound.bestAlign, 1))
                std::cout << c;
            std::cout << std::endl;
            //--------------------
        }

        if (checkEraseV)
        {
            for (int i = eraseV.size() - 1; i >= 0; --i)
            {
                if(eraseV[i])
                {
                    goldRnaAligns.push_back(rnaAligns[i]);
                    rnaAligns.erase(rnaAligns.begin() + i);
                    erase(alignsSimd, i);
                    erase(resultsSimd, i);
                    eraseV.erase(eraseV.begin() + i);
                }
            }
        }
        //if (options.verbose > 0) std::cerr << "|";

    }
    if (options.verbose > 0) std::cerr << std::endl;
    _VV(options, "map computation time = " << boutime << "\nlemon MWM time       = " << lemtime
                                           << "\ngreedy MWM time      = " << mwmtime);

// timer stop
    std::chrono::steady_clock::time_point endChrono= std::chrono::steady_clock::now();
    std::clock_t end = std::clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

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

    rnaAligns.insert( rnaAligns.end(), goldRnaAligns.begin(), goldRnaAligns.end() );

    if(rnaAligns.size() > 0) // This is a multiple alignment an the T-Coffee library must be printed
    {
        createTCoffeeLib(options, filecontents1, filecontents2, rnaAligns);
    }




// Print elapsed time
    _VV(options, "\nTime difference chrono = " << std::chrono::duration_cast<std::chrono::seconds>(endChrono - beginChrono).count()); //std::chrono::microseconds
    _VV(options, "\nTime difference = " << elapsed_secs);

    return 0;
}

