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

void createMask(TRnaAlign & rnaAlign, TAlign const & align)
{
    Row<TAlign>::Type row0 = row(align, 0);
    Row<TAlign>::Type row1 = row(align, 1);

    clear(rnaAlign.mask);

    for (std::size_t column = 0u; column < length(row0); ++column)
    {
        if (!isGap(row0, column) && !isGap(row1, column))
            appendValue(rnaAlign.mask, std::make_pair(toSourcePosition(row0, column), toSourcePosition(row1, column)));
    }
}

// ----------------------------------------------------------------------------
// Function computeUpperBound()
// ----------------------------------------------------------------------------

void computeUpperBoundScore(TRnaAlign & rnaAlign)
{
    TScoreValue upperBoundScore = 0;
    rnaAlign.slm = 0;
    for (unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        if (rnaAlign.upperBoundVect[i].maxProbScoreLine > 0)
        {
            upperBoundScore += rnaAlign.upperBoundVect[i].maxProbScoreLine;
            std::cerr << "  " << rnaAlign.upperBoundVect[i].maxProbScoreLine;
            if (rnaAlign.upperBoundVect[i].seq1IndexPairLine !=
                rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq2IndexPairLine].seq1Index)
            {
                // the edges are not paired
                ++rnaAlign.slm;
                std::cerr << "*";
            }
        }
    }
    std::cerr << std::endl;
    rnaAlign.upperBound = upperBoundScore;
};

// ----------------------------------------------------------------------------
// Function computeBound()
// ----------------------------------------------------------------------------

void computeLowerAndUpperBoundScore(TRnaAlign & rnaAlign)
{
    TScoreValue sumU = 0;
    TScoreValue sumL = 0;
    rnaAlign.slm = 0;
    for(unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i) {
        if (rnaAlign.upperBoundVect[i].maxProbScoreLine > 0) {
            // the edges are not paired
            if (rnaAlign.upperBoundVect[i].seq1Index !=
                rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq2IndexPairLine].seq1IndexPairLine)
            {
                sumU += rnaAlign.upperBoundVect[i].maxProbScoreLine;
                ++rnaAlign.slm;
            }
            else
            {
                sumL += rnaAlign.upperBoundVect[i].maxProbScoreLine;
            }
//            std::cout << "Pairs " << rnaAlign.upperBoundVect[i].seq1Index << ":";
//            std::cout << rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq2IndexPairLine].seq1IndexPairLine << std::endl;
//            std::cout << rnaAlign.upperBoundVect[i].maxProbScoreLine << "\t";
//            std::cout << rnaAlign.upperBoundVect[i].seq1Index << ":" << i << "\t";
//            std::cout << rnaAlign.upperBoundVect[i].seq1IndexPairLine << ":";
//            std::cout << rnaAlign.upperBoundVect[i].seq2IndexPairLine << std::endl;
        }
    }
    rnaAlign.upperBound = sumU + sumL;
    rnaAlign.lowerBound = sumL;
//    std::cout << "upperBound = " << sumU + sumL << std::endl;
//    std::cout << "lowerBound = " << sumL << std::endl;
};

// ----------------------------------------------------------------------------
// Function computeBounds() version that make use of the lemon MWM
// ----------------------------------------------------------------------------

void computeBounds(TRnaAlign & rnaAlign, TMapVect * lowerBound4Lemon) // upper bound computation
{
    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;

    // Clear the maxProbScoreLine of the upper bound
    for (std::size_t idx = 0; idx < length(rnaAlign.upperBoundVect); ++idx)
    {
        rnaAlign.upperBoundVect[idx].maxProbScoreLine = 0; // reset best line score
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
                if (rnaAlign.upperBoundVect[line.second].maxProbScoreLine < edgeWeight / 2.0)
                {
                    rnaAlign.upperBoundVect[line.second].maxProbScoreLine = edgeWeight / 2.0;
                    rnaAlign.upperBoundVect[line.second].seq1Index = line.first;
                    rnaAlign.upperBoundVect[line.second].seq1IndexPairLine = value(adj_it1);
                    rnaAlign.upperBoundVect[line.second].seq2IndexPairLine = value(adj_it2);
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
        rnaAlign.forMinBound.upperBoundVect = rnaAlign.upperBoundVect;
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
        rnaAlign.forScore.upperBoundVect = rnaAlign.upperBoundVect;
        rnaAlign.forScore.mask = rnaAlign.mask;
    }
}

void saveBestAligns(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
    saveBestAlignMinBound(rnaAlign, align, alignScore, index);
    saveBestAlignScore(rnaAlign, align, alignScore, index);
}


void updateLambda(TRnaAlign & rnaAlign) {
    for (size_t i = 0; i < length(rnaAlign.upperBoundVect); ++i) {
        struct boundStruct const & ub = rnaAlign.upperBoundVect[i];

        if (ub.maxProbScoreLine > 0) { // skip if no interactions
            // get lambda struct of paired seq1 position
            struct lambWeightStruct & lambWeight = rnaAlign.lamb[ub.seq1Index].map[i];

            // the interaction edges are not paired
            if (ub.seq1Index != rnaAlign.upperBoundVect[ub.seq2IndexPairLine].seq1IndexPairLine)
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
