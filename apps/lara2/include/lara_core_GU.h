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
// Function fillBound()
// ----------------------------------------------------------------------------

void fillBound(TBound & BoundVect, TRnaAlign const & rnaAlign, unsigned x, unsigned y, TScoreValue halfEdgesProb)
{
    BoundVect[rnaAlign.mask[x].second].maxProbScoreLine = halfEdgesProb;
    BoundVect[rnaAlign.mask[x].second].seq1Index = rnaAlign.mask[x].first;
    BoundVect[rnaAlign.mask[x].second].seq1IndexPairLine = rnaAlign.mask[y].first;
    BoundVect[rnaAlign.mask[x].second].seq2IndexPairLine = rnaAlign.mask[y].second;
};

// ----------------------------------------------------------------------------
// Function computeUpperBound()
// ----------------------------------------------------------------------------

void computeUpperBoundScore(TRnaAlign & rnaAlign)
{
    TScoreValue sum = 0;
    rnaAlign.slm = 0;
    for (unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        if (rnaAlign.upperBoundVect[i].maxProbScoreLine > 0)
        {
            sum += rnaAlign.upperBoundVect[i].maxProbScoreLine;
            if (rnaAlign.upperBoundVect[i].seq1Index !=
                rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq2IndexPairLine].seq1IndexPairLine)
            {
                // the edges are not paired
                ++rnaAlign.slm;
            }
        }
    }
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
    TScoreValue edgesProb, halfEdgesProb;

    // Clear the maxProbScoreLine of the upper bound
    for (std::size_t idx = 0; idx < length(rnaAlign.upperBoundVect); ++idx)
    {
        rnaAlign.upperBoundVect[idx].maxProbScoreLine = 0; // reset best line score
    }

    for (unsigned i = 0; i < length(rnaAlign.mask) - 1; ++i)
    {
        std::pair<unsigned, unsigned> line = rnaAlign.mask[i];
        //String<unsigned> adjVect1;
        // get all interacting positions for current mask[i]
        //getVertexAdjacencyVector(adjVect1, graph1, rnaAlign.mask[i].first);
        //for (unsigned j = 0; j < length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)

        for (RnaAdjacencyIterator adj_it1(graph1, line.first); !atEnd(adj_it1); goNext(adj_it1))
        {
            double edgeWeight1 = cargo(findEdge(graph1, line.first, value(adj_it1)));
            for (RnaAdjacencyIterator adj_it2(graph2, line.second); !atEnd(adj_it2); goNext(adj_it2))
            {
                double edgeWeight = edgeWeight1 + cargo(findEdge(graph2, line.second, value(adj_it2)));
                if (lowerBound4Lemon != NULL)
                    (*lowerBound4Lemon)[value(adj_it1)][value(adj_it2)] = edgeWeight;

                std::cerr << "Interaction match: " << line.first << " - " << line.second
                          << " | " << value(adj_it1) << " - " << value(adj_it2) << "\tprob = "
                          << edgeWeight1 << " + " << edgeWeight-edgeWeight1 << "\n";

                if (rnaAlign.upperBoundVect[line.second].maxProbScoreLine < edgeWeight / 2.0)
                {
                    //fillBound(rnaAlign.upperBoundVect, rnaAlign, i, w, halfEdgesProb);
                    rnaAlign.upperBoundVect[line.second].maxProbScoreLine = edgeWeight / 2.0;
                    rnaAlign.upperBoundVect[line.second].seq1Index = line.first;
                    rnaAlign.upperBoundVect[line.second].seq1IndexPairLine = value(adj_it1);
                    rnaAlign.upperBoundVect[line.second].seq2IndexPairLine = value(adj_it2);
                    std::cerr << "updated\n";
                }
            }

/*
            if (value(adj_it1) < rnaAlign.mask[i].first) //TODO Are values of adj_it sorted? Then move cond. to loop
                continue; // detect interactions only to the right

            for (unsigned w = i + 1; w < length(rnaAlign.mask); ++w)
            {
                if (value(adj_it1) != rnaAlign.mask[w].first)
                    continue; // look for mask element mask[w] that is adjacent
                // TODO here we loose important edges!

                //String<unsigned> adjVect2;
                //getVertexAdjacencyVector(adjVect2, graph2, rnaAlign.mask[i].second);
                //for (unsigned z = 0; z < length(adjVect2) && adjVect2[z] >= rnaAlign.mask[i].second; ++z)
                for (RnaAdjacencyIterator adj_it2(graph2, rnaAlign.mask[i].second); !atEnd(adj_it2); goNext(adj_it2))
                {
                    if (value(adj_it2) < rnaAlign.mask[i].second)
                        continue; // detect interactions only to the right

                    if (value(adj_it2) == rnaAlign.mask[w].second) // found interaction match!
                    {
                        std::cerr << "Interaction match: " << rnaAlign.mask[i].first << " - " << rnaAlign.mask[w].first
                                  << " | " << rnaAlign.mask[i].second << " - " << rnaAlign.mask[w].second << "\tprob = "
                                  << cargo(findEdge(graph1, rnaAlign.mask[i].first, rnaAlign.mask[w].first)) << " + "
                                  << cargo(findEdge(graph2, rnaAlign.mask[i].second, rnaAlign.mask[w].second)) << "\n";
                        edgesProb = (cargo(findEdge(graph1, rnaAlign.mask[i].first, rnaAlign.mask[w].first)) +
                                     cargo(findEdge(graph2, rnaAlign.mask[i].second, rnaAlign.mask[w].second)));

                        // Add edge for the LowerBoundGraph
                        if (lowerBound4Lemon != NULL)
                            (*lowerBound4Lemon)[i][w] = edgesProb;

                        // Compute the half edge probability to be used for the upper bound level
                        halfEdgesProb = edgesProb / 2.0;
                        if (rnaAlign.upperBoundVect[rnaAlign.mask[i].second].maxProbScoreLine < halfEdgesProb)
                        {
                            fillBound(rnaAlign.upperBoundVect, rnaAlign, i, w, halfEdgesProb);
                        }
                        if (rnaAlign.upperBoundVect[rnaAlign.mask[w].second].maxProbScoreLine < halfEdgesProb )
                        {
                            fillBound(rnaAlign.upperBoundVect, rnaAlign, w, i, halfEdgesProb);
                        }
                    }
                }
            }
            */

        }
    }
    // computeUpperBound(rnaAlign);
}

// ----------------------------------------------------------------------------
// Function computeBounds() version that make use of an approximation of the MWM
// ----------------------------------------------------------------------------

/*
//template <typename TRnaAlign>
void computeBounds(TRnaAlign & rnaAlign)
{
    RnaStructureGraph & graph1 = rnaAlign.bppGraphH;
    RnaStructureGraph & graph2 = rnaAlign.bppGraphV;
    TScoreValue edgesProb, halfEdgesProb;
//  Clear the maxProbScoreLine of the upper bound
    for (unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        rnaAlign.upperBoundVect[i].maxProbScoreLine = 0;
    }
    for (unsigned i = 0; i < length(rnaAlign.mask) - 1; ++i)
    {
        String<unsigned> adjVect1;
        getVertexAdjacencyVector(adjVect1, graph1.inter, rnaAlign.mask[i].first);
        for (unsigned j = 0; j < length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)
        {
            for (unsigned w = i + 1; w < length(rnaAlign.mask); ++w)
            {
                if (adjVect1[j] == rnaAlign.mask[w].first)
                {
                    String<unsigned> adjVect2;
                    getVertexAdjacencyVector(adjVect2, graph2.inter, rnaAlign.mask[i].second);
                    for (unsigned z = 0; z < length(adjVect2) && adjVect2[z] >= rnaAlign.mask[i].second; ++z)
                    {
                        if (adjVect2[z] == rnaAlign.mask[w].second)
                        {
                            edgesProb = (cargo(findEdge(graph1.inter, rnaAlign.mask[i].first,
                                                        rnaAlign.mask[w].first)) +
                                         cargo(findEdge(graph2.inter, rnaAlign.mask[i].second,
                                                        rnaAlign.mask[w].second)));
                            // Compute the half edge probability to be used for the upper bound level
                            halfEdgesProb = edgesProb / 2.0;
                            if (rnaAlign.upperBoundVect[rnaAlign.mask[i].second].maxProbScoreLine < halfEdgesProb)
                            {
                                fillBound(rnaAlign.upperBoundVect, rnaAlign, i, w, halfEdgesProb);
                            }
                            if (rnaAlign.upperBoundVect[rnaAlign.mask[w].second].maxProbScoreLine < halfEdgesProb )
                            {
                                fillBound(rnaAlign.upperBoundVect, rnaAlign, w, i, halfEdgesProb);
                            }
                        }
                    }
                }
            }
        }
    }
    // computeBound(rnaAlign);
}


// ----------------------------------------------------------------------------
// Function computeBoundsTest()
// ----------------------------------------------------------------------------

//template <typename TOption, typename TRnaAlign, typename TMapVect>
void computeBoundsTest(TRnaAlign & rnaAlign, TMapVect & lowerBound4Lemon)
{
    RnaStructureGraph & graph1 = rnaAlign.bppGraphH;
    RnaStructureGraph & graph2 = rnaAlign.bppGraphV;
    TScoreValue edgesProb, halfEdgesProb;
//  Clear the maxProbScoreLine of the upper bound
    for(unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        rnaAlign.upperBoundVect[i].maxProbScoreLine = 0;
//        std::cout << rnaAlign.upperBoundVect[i].maxProbScoreLine << "\t";
    }
    std::cout << std::endl;
    clearVertices(rnaAlign.lowerBoundGraph);
    // Add vertex for the LowerBoundGraph
    for(unsigned i = 0; i < length(rnaAlign.mask); ++i)
    {
        addVertex(rnaAlign.lowerBoundGraph);
    }
    unsigned ll = 0;
    for (unsigned i = 0; i < length(rnaAlign.mask) - 1; ++i)
    {
        String<unsigned> adjVect1;
        getVertexAdjacencyVector(adjVect1, graph1.inter, rnaAlign.mask[i].first);
        for(unsigned j = 0; j< length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)
        {
            for (unsigned w = i + 1; w < length(rnaAlign.mask); ++w)
            {
                if (adjVect1[j] == rnaAlign.mask[w].first)
                {
                    std::cout << "Couple of nt at position " << i << " of mask = " << rnaAlign.mask[i].first << ":"
                              << rnaAlign.mask[i].second << std::endl;
                    std::cout << "Couple of nt at position " << w << " of mask = " << rnaAlign.mask[w].first << ":"
                              << rnaAlign.mask[w].second << std::endl;
                    std::cout << "Probability on seq1 " << rnaAlign.mask[i].first << ":" << rnaAlign.mask[w].first
                              << " = " << cargo(findEdge(graph1.inter, rnaAlign.mask[i].first,
                                                         rnaAlign.mask[w].first)) << std::endl;
                    String<unsigned> adjVect2;
                    getVertexAdjacencyVector(adjVect2, graph2.inter, rnaAlign.mask[i].second);
                    for (unsigned z = 0; z < length(adjVect2) && adjVect2[z] >= rnaAlign.mask[i].second; ++z)
                    {
                        if (adjVect2[z] == rnaAlign.mask[w].second)
                        {
                            std::cout << "Probability on seq2 " << rnaAlign.mask[w].second << ":"
                                      << rnaAlign.mask[i].second << " = "
                                      << cargo(findEdge(graph2.inter, rnaAlign.mask[i].second,
                                                        rnaAlign.mask[w].second)) << std::endl;
                            edgesProb = (cargo(findEdge(graph1.inter, rnaAlign.mask[i].first, rnaAlign.mask[w].first)) +
                                         cargo(findEdge(graph2.inter, rnaAlign.mask[i].second, rnaAlign.mask[w].second)));
                            // Add edge for the LowerBoundGraph
                            addEdge(rnaAlign.lowerBoundGraph, i, w, edgesProb);
                            lowerBound4Lemon[i][w] = edgesProb;
                            ++ll;
                            // Compute the half edge probability to be used for the upper bound level
                            halfEdgesProb = edgesProb / 2.0;
                            // std::cout << "Sum of probabilities = " << edgesProb << std::endl;

                            if (rnaAlign.upperBoundVect[rnaAlign.mask[i].second].maxProbScoreLine < halfEdgesProb)
                            {
                                std::cout << "Previous prob of pair " << rnaAlign.mask[i].second << " = "
                                          << rnaAlign.upperBoundVect[rnaAlign.mask[i].second].maxProbScoreLine
                                          << std::endl;
                                fillBound(rnaAlign.upperBoundVect, rnaAlign, i, w, halfEdgesProb);
                                std::cout << "Updated prob of pair " << rnaAlign.mask[i].second << " = "
                                          << rnaAlign.upperBoundVect[rnaAlign.mask[i].second].maxProbScoreLine
                                          << std::endl;
                            }
                            if ( rnaAlign.upperBoundVect[rnaAlign.mask[w].second].maxProbScoreLine < halfEdgesProb )
                            {
                                std::cout << "Previous prob of pair " << rnaAlign.mask[w].second << " = "
                                          << rnaAlign.upperBoundVect[rnaAlign.mask[w].second].maxProbScoreLine
                                          << std::endl;
                                fillBound(rnaAlign.upperBoundVect, rnaAlign, w, i, halfEdgesProb);
                                std::cout << "Updated prob of pair " << rnaAlign.mask[w].second << " = "
                                          << rnaAlign.upperBoundVect[rnaAlign.mask[w].second].maxProbScoreLine
                                          << std::endl;
                            }
                            std::cout << adjVect2[z] << "\t";
                            std::cout << rnaAlign.mask[i].first << ":" << rnaAlign.mask[i].second << std::endl;
                            std::cout << rnaAlign.mask[w].first << ":" << rnaAlign.mask[w].second << std::endl;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Number of edges = " << ll << std::endl;
    std::cout << "Upper bound vector of size " << length(rnaAlign.upperBoundVect) << std::endl;

    computeBound(rnaAlign);

    std::cout << "upperBound = " << rnaAlign.upperBound << std::endl;
    std::cout << "lowerBound = " << rnaAlign.lowerBound << std::endl;
    std::cout << "lowerBound from lemon::MWM primal = " << rnaAlign.lowerLemonBound.mwmPrimal << " dual = "
              << rnaAlign.lowerLemonBound.mwmDual << std::endl;
    std::cout << rnaAlign.lowerBoundGraph << std::endl;
}
*/

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
