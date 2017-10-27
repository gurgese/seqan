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
// This file contains the variable definition and structures of laragu
// application.
// ==========================================================================

#ifndef _INCLUDE_ALIGNMENT_EDGES_H_
#define _INCLUDE_ALIGNMENT_EDGES_H_

#include <map>
#include <utility>
#include <seqan/graph_algorithms.h>
#include "lemon_graph.h"

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function createMask()
// ----------------------------------------------------------------------------

template <typename TOptions>
void createMask(TRnaAlign & rnaAlign, TAlign const & align, TOptions const & options)
{
    typedef typename Iterator<Gaps<TSequence, seqan::ArrayGaps> const, Standard>::Type TGapsIter;

    Row<TAlign>::Type row0 = row(align, 0);
    Row<TAlign>::Type row1 = row(align, 1);

    clear(rnaAlign.mask);
    rnaAlign.sequenceScore = 0;

    // Get iterators.
    TGapsIter it0      = begin(row0);
    TGapsIter itEnd0   = end(row0);
    TGapsIter it1      = begin(row1);
    TGapsIter itEnd1   = end(row1);

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;

    for (unsigned column = 0u; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1, ++column)
    {
        // Gaps in first sequence
        if (isGap(it0))
        {
            if (!isGapOpen0)
            {
                rnaAlign.sequenceScore += options.laraGapOpen;
            }
            else
            {
                rnaAlign.sequenceScore += options.laraGapExtend;
            }
            isGapOpen0 = true;
        }
        else
        {
            isGapOpen0 = false;
        }

        // Gaps in second sequence
        if (isGap(it1))
        {
            if (!isGapOpen1)
            {
                rnaAlign.sequenceScore += options.laraGapOpen;
            }
            else
            {
                rnaAlign.sequenceScore += options.laraGapExtend;
            }
            isGapOpen1 = true;
        }
        else
        {
            isGapOpen1 = false;
        }

        // Match or mismatch
        if (!isGap(it0) && !isGap(it1))
        {
            rnaAlign.sequenceScore += score(options.laraScoreMatrix, *it0, *it1);

            // create mask entry
            appendValue(rnaAlign.mask, std::make_pair(toSourcePosition(row0, column), toSourcePosition(row1, column)));
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);
}

// ----------------------------------------------------------------------------
// Function computeUpperBound()
// ----------------------------------------------------------------------------

void computeUpperBoundScore(TRnaAlign & rnaAlign)
{
    TScoreValue sum = 0;
    rnaAlign.slm = 0;
    for (unsigned i = 0; i < length(rnaAlign.weightLineVect); ++i)
    {
        if (rnaAlign.weightLineVect[i].maxProbScoreLine > 0)
        {
            sum += rnaAlign.weightLineVect[i].maxProbScoreLine;
            std::cerr << "  " << rnaAlign.weightLineVect[i].maxProbScoreLine;
            if (rnaAlign.weightLineVect[i].seq1Index !=
                rnaAlign.weightLineVect[rnaAlign.weightLineVect[i].seq2IndexPairLine].seq1IndexPairLine)
            {
                // the edges are not paired
                ++rnaAlign.slm;
                std::cerr << "*";
            }
        }
    }
    std::cerr << std::endl;
    rnaAlign.upperBound = sum;
};

// ----------------------------------------------------------------------------
// Function computeBound()
// ----------------------------------------------------------------------------

void computeLowerAndUpperBoundScore(TRnaAlign & rnaAlign)
{
    TScoreValue sumU = 0;
    TScoreValue sumL = 0;
    rnaAlign.slm = 0;
    for(unsigned i = 0; i < length(rnaAlign.weightLineVect); ++i) {
        if (rnaAlign.weightLineVect[i].maxProbScoreLine > 0) {
            // the edges are not paired
            if (rnaAlign.weightLineVect[i].seq1Index !=
                rnaAlign.weightLineVect[rnaAlign.weightLineVect[i].seq2IndexPairLine].seq1IndexPairLine)
            {
                sumU += rnaAlign.weightLineVect[i].maxProbScoreLine;
                ++rnaAlign.slm;
            }
            else
            {
                sumL += rnaAlign.weightLineVect[i].maxProbScoreLine;
            }
//            std::cout << "Pairs " << rnaAlign.weightLineVect[i].seq1Index << ":";
//            std::cout << rnaAlign.weightLineVect[rnaAlign.weightLineVect[i].seq2IndexPairLine].seq1IndexPairLine << std::endl;
//            std::cout << rnaAlign.weightLineVect[i].maxProbScoreLine << "\t";
//            std::cout << rnaAlign.weightLineVect[i].seq1Index << ":" << i << "\t";
//            std::cout << rnaAlign.weightLineVect[i].seq1IndexPairLine << ":";
//            std::cout << rnaAlign.weightLineVect[i].seq2IndexPairLine << std::endl;
        }
    }
    rnaAlign.upperBound = sumU + sumL;
    rnaAlign.lowerBound = sumL;
//    std::cout << "upperBound = " << sumU + sumL << std::endl;
//    std::cout << "lowerBound = " << sumL << std::endl;
};


void computeLowerBound(TRnaAlign & rnaAlign, TMapVect * lowerBound4Lemon)
{
//    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
//    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;
    // iterate all lines that are present in the alignment
    for (std::pair<unsigned, unsigned> const & line : rnaAlign.mask)
    {
        if (rnaAlign.lamb[line.first].map.count(line.second) > 0)
        {
            lambWeightStruct &lamb = rnaAlign.lamb[line.first].map[line.second];
            std::cout << lamb.seq1IndexPairLine << std::endl;
            if (lamb.seq1IndexPairLine == rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].seq1IndexInter &&
                lamb.seq2IndexPairLine == rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].seq2IndexInter &&
                !rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].fromUBPairing) {
//                std::cout << lamb.seq1IndexPairLine << " == "
//                          << rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].seq1IndexInter << " && " <<
//                          lamb.seq2IndexPairLine << " == "
//                          << rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].seq2IndexInter << " &&" <<
//                          !rnaAlign.lamb[lamb.seq1IndexInter].map[lamb.seq2IndexInter].fromUBPairing << std::endl;
                (*lowerBound4Lemon)[lamb.seq1IndexPairLine][lamb.seq1IndexInter] = lamb.maxProbScoreLine * 2;
            } else
            {
                ++rnaAlign.slm;
            }
        }
    }
    std::cerr << "The sequence score is " << rnaAlign.sequenceScore << std::endl;
//    std::cout << rnaAlign.slm << std::endl;
}

// ----------------------------------------------------------------------------
// Function computeBounds() version that make use of the lemon MWM
// ----------------------------------------------------------------------------

void computeBounds(TRnaAlign & rnaAlign, TMapVect * lowerBound4Lemon) // upper bound computation
{
    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;

    // Clear the maxProbScoreLine of the upper bound
    for (std::size_t idx = 0; idx < length(rnaAlign.weightLineVect); ++idx)
    {
        rnaAlign.weightLineVect[idx].maxProbScoreLine = 0; // reset best line score
    }

    // iterate all lines that are present in the alignment
    for (std::pair<unsigned, unsigned> const & line : rnaAlign.mask)
    {
        // outgoing interactions in the first sequence
        for (RnaAdjacencyIterator adj_it1(graph1, line.first); !atEnd(adj_it1); goNext(adj_it1))
        {
            double edgeWeight1 = cargo(findEdge(graph1, line.first, value(adj_it1)));

            // outgoing interactions in the second sequence
            for (RnaAdjacencyIterator adj_it2(graph2, line.second); !atEnd(adj_it2); goNext(adj_it2))
            {
                double edgeWeight = edgeWeight1 + cargo(findEdge(graph2, line.second, value(adj_it2)));

                // for lower bound the directed interactions must occur in both orientations
                if (lowerBound4Lemon != NULL && line.first < value(adj_it1))
                {
                    for (std::pair<unsigned, unsigned> const & pairline : rnaAlign.mask)
                    { //TODO make this loop more efficient
                        if (pairline.first == value(adj_it1) && pairline.second == value(adj_it2))
                            (*lowerBound4Lemon)[line.first][pairline.first] = edgeWeight;
                    }
                }

                std::cerr << "Interaction match: " << line.first+1 << " - " << line.second+1
                          << " | " << value(adj_it1)+1 << " - " << value(adj_it2)+1 << "\tprob = "
                          << edgeWeight1 << " + " << edgeWeight-edgeWeight1 << "\n";

                // for upper bound do not care if interactions are closed by a line
                if (rnaAlign.weightLineVect[line.second].maxProbScoreLine < edgeWeight / 2.0)
                {
                    rnaAlign.weightLineVect[line.second].maxProbScoreLine = edgeWeight / 2.0;
                    rnaAlign.weightLineVect[line.second].seq1Index = line.first;
                    rnaAlign.weightLineVect[line.second].seq1IndexPairLine = value(adj_it1);
                    rnaAlign.weightLineVect[line.second].seq2IndexPairLine = value(adj_it2);
                    std::cerr << "updated\n";
                }
            }
        }
    }
}

void computeLowerBoundGreedy(TMapVect & interactions, TRnaAlign & rnaAlign)
{
    TLowerBoundGraph graph;

    // add vertices
    forEach(interactions, [&graph] (TMap const &) { addVertex(graph); });

    // add edges
    for (unsigned vertexIdx = 0; vertexIdx < length(interactions); ++vertexIdx)
        for (auto edgeCargo = interactions[vertexIdx].begin(); edgeCargo != interactions[vertexIdx].end(); ++edgeCargo)
            addEdge(graph, vertexIdx, edgeCargo->first, edgeCargo->second);

    rnaAlign.lowerGreedyBound = maximumWeightedMatchingGreedy<5>(graph);
};

// ----------------------------------------------------------------------------
// Function saveBestAlign()
// ----------------------------------------------------------------------------

void saveBestAlignMinBound(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, unsigned index)
{
//    if ((rnaAlign.upperBound - rnaAlign.lowerBound) < (rnaAlign.upperMinBound - rnaAlign.lowerMinBound))
    if (rnaAlign.stepSize < rnaAlign.forMinBound.stepSizeBound)  //TODO check if this <= is expensive
    {
        std::cerr << "update best min bound" << std::endl;
        rnaAlign.forMinBound.it = index; //to be used for the best lower bound
        rnaAlign.forMinBound.lowerBound = rnaAlign.lowerBound;
        rnaAlign.forMinBound.upperBound = rnaAlign.upperBound;
        rnaAlign.forMinBound.stepSizeBound = rnaAlign.stepSize;
        rnaAlign.forMinBound.bestAlign = align;
        rnaAlign.forMinBound.bestAlignScore = alignScore;
        rnaAlign.forMinBound.weightLineVect = rnaAlign.weightLineVect;
        rnaAlign.forMinBound.mask = rnaAlign.mask;
    }
}
void saveBestAlignScore(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
//    if ((rnaAlign.upperBound - rnaAlign.lowerBound) < (rnaAlign.upperMinBound - rnaAlign.lowerMinBound))
    if (rnaAlign.forScore.bestAlignScore < alignScore)
    {
        std::cerr << "update best score" << std::endl;
        rnaAlign.forScore.it = index; //to be used for the best lower bound
        rnaAlign.forScore.lowerBound = rnaAlign.lowerBound;
        rnaAlign.forScore.upperBound = rnaAlign.upperBound;
        rnaAlign.forScore.stepSizeBound = rnaAlign.stepSize;
        rnaAlign.forScore.bestAlign = align;
        rnaAlign.forScore.bestAlignScore = alignScore;
        rnaAlign.forScore.weightLineVect = rnaAlign.weightLineVect;
        rnaAlign.forScore.mask = rnaAlign.mask;
    }
}

void saveBestAligns(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
    saveBestAlignMinBound(rnaAlign, align, alignScore, index);
    saveBestAlignScore(rnaAlign, align, alignScore, index);
}

template <typename TValueScoreLine>
void updateLambdaLine(TValueScoreLine & maxProbScoreLine, unsigned & seqIndexInter, Graph<Undirected<double> > const & graph, unsigned const & position)
{
    for (RnaAdjacencyIterator adj_it(graph, position); !atEnd(adj_it); goNext(adj_it))
    {
        double edgeWeight = cargo(findEdge(graph, position, value(adj_it)));
//        std::cout << line.first << " : " << value(adj_it1) << " \ " << edgeWeight << " - " << line.second << std::endl;
        if(maxProbScoreLine < edgeWeight)
        {
            maxProbScoreLine = edgeWeight;
            seqIndexInter = value(adj_it);
        }
    }
}

void switchIndex (unsigned & index1, unsigned & index2)
{
    unsigned tmp = index1;
    index1 = index2;
    index2 = tmp;
}


void computeLambda(TRnaAlign & rnaAlign, bool const & saveFoundInterPair)
{
    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;

    // iterate all lines that are present in the alignment
    for (std::pair<unsigned, unsigned> const & line : rnaAlign.mask)
    {
        std::cout << degree(graph1, line.first) << " - " << degree(graph2, line.second) << std::endl;
        // outgoing interactions in the first sequence
        bool proceed = false;
        if (degree(graph1, line.first) > 0 && degree(graph2, line.second) > 0)
        {
            if (rnaAlign.lamb[line.first].map.count(line.second) == 0)
            {
                proceed = true;
            } else if (rnaAlign.lamb[line.first].map[line.second].fromUBPairing)
            {
                proceed = true;
                rnaAlign.lamb[line.first].map[line.second].fromUBPairing = false;
            }
            if (proceed)
            {
                lambWeightStruct & lambda = rnaAlign.lamb[line.first].map[line.second];
                updateLambdaLine(lambda.maxProbScoreLine1, lambda.seq1IndexInter, graph1, line.first);
                lambda.seq1IndexPairLine = line.first;
//            std::cout << lambda.seq1IndexPairLine << " : " << lambda.seq1IndexInter << " | " << lambda.maxProbScoreLine1
//                      << " - " << line.second << std::endl;
                updateLambdaLine(lambda.maxProbScoreLine2, lambda.seq2IndexInter, graph2, line.second);
                lambda.seq2IndexPairLine = line.second;
//            std::cout << lambda.seq2IndexPairLine << " : " << lambda.seq2IndexInter << " | " << lambda.maxProbScoreLine2
//                      << " - " << line.first << std::endl;
                lambda.maxProbScoreLine = (lambda.maxProbScoreLine1 + lambda.maxProbScoreLine2) / 2;

                std::cout << lambda.seq1IndexPairLine << " : " << lambda.seq2IndexPairLine << " = "
                          << lambda.maxProbScoreLine << " | " << lambda.seq1IndexInter << " : "
                          << lambda.seq2IndexInter << " (" << lambda.fromUBPairing << ")" << std::endl;
                rnaAlign.lamb[line.first].map[line.second].fromUBPairing = false;
                //            std::cout << rnaAlign.lamb[lambda.seq1IndexInter].map.count(lambda.seq1IndexPairLine) << std::endl;
                if (saveFoundInterPair && rnaAlign.lamb[lambda.seq1IndexInter].map.count(lambda.seq1IndexPairLine) == 0)
                    // if the map is not empty means that the line have been already evaluated in a previous iteration
                {
                    lambWeightStruct & lambda = rnaAlign.lamb[line.first].map[line.second];
                    struct lambWeightStruct &lambdaPair = rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq1IndexPairLine];
                    lambdaPair = rnaAlign.lamb[line.first].map[line.second];
                    switchIndex(lambdaPair.seq1IndexPairLine, lambdaPair.seq1IndexInter);
                    switchIndex(lambdaPair.seq2IndexPairLine, lambdaPair.seq2IndexInter);
                    lambdaPair.fromUBPairing = true;

                    std::cout << lambdaPair.seq1IndexPairLine << " : " << lambdaPair.seq2IndexPairLine << " = "
                              << lambdaPair.maxProbScoreLine << " | " << lambdaPair.seq1IndexInter << " : "
                              << lambdaPair.seq2IndexInter << " (" << lambdaPair.fromUBPairing << ")" << std::endl;
                }
            }
        }
//        rnaAlign[line.first].map[line.second].maxProbScoreLine = lambWeight/2;
    }
}

void updateLambda(TRnaAlign & rnaAlign)  // TODO REMOVE THIS FUNCTION
{
    for (size_t i = 0; i < length(rnaAlign.weightLineVect); ++i)
    {
        struct weightLineStruct const & ub = rnaAlign.weightLineVect[i];

        if (ub.maxProbScoreLine > 0) { // skip if no interactions
            // get lambda struct of paired seq1 position
            struct lambWeightStruct & lambWeight = rnaAlign.lamb[ub.seq1Index].map[i];

            // the interaction edges are not paired
            if (ub.seq1Index != rnaAlign.weightLineVect[ub.seq2IndexPairLine].seq1IndexPairLine)
            {
                // if C < D add stepSize to lambda
                if (ub.seq1Index < ub.seq1IndexPairLine)
                {
// TODO check if this strategy is properly working a positive score is assigned to the left-side alignments.
// Maybe a double side strategy should be tested
                    lambWeight.step += rnaAlign.stepSize;
                    // Note, the default initializer is callet the fist time that set the value to 0
                } else {
                    lambWeight.step -= rnaAlign.stepSize;
                    // Note, the default initializer is callet the fist time that set the value to 0
                }
            }
// Save the maximum interaction weight to be used for the computation of profit of a line
            if (lambWeight.maxProbScoreLine < ub.maxProbScoreLine)
            {
                std::cerr << "lambda[" << i << "] changes: " << lambWeight.seq1IndexPairLine << ","
                          << lambWeight.seq2IndexPairLine << "->" << ub.seq1IndexPairLine << ","
                          << ub.seq2IndexPairLine << " \tnew value = " << ub.maxProbScoreLine << std::endl;
                lambWeight.maxProbScoreLine = ub.maxProbScoreLine;
                lambWeight.seq1IndexPairLine = ub.seq1IndexPairLine;
                lambWeight.seq2IndexPairLine = ub.seq2IndexPairLine;
            }
        }
    }
}
#endif //_INCLUDE_ALIGNMENT_EDGES_H_
