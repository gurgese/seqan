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

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function saveBestAlign()
// ----------------------------------------------------------------------------

//template <typename TOption, typename TAlign, typename TScoreValue, typename TRnaAlign>
void saveBestAlign(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore)
{
    if (rnaAlign.bestAlignScore < alignScore)
    {
        rnaAlign.bestAlign = align;
        rnaAlign.bestAlignScore = alignScore;
    }
};

// ----------------------------------------------------------------------------
// Function maskCreator()
// ----------------------------------------------------------------------------

//template <typename TOption, typename TAlign, typename TScoreValue, typename TRnaAlign>
void maskCreator(TRnaAlign & rnaAlign, TAlign const & align)
{
    unsigned row0Begin = clippedBeginPosition(row(align, 0));
    unsigned row1Begin = clippedBeginPosition(row(align, 1));
    unsigned gap0 = 0;
    unsigned gap1 = 0;
    unsigned j = 0;
//    std::cout << "row0Begin " << row0Begin << " row1Begin " << row1Begin << std::endl;
    //TODO this function can be simplifyed using the function toViewPosition(row0, i)
    // Initialize the mask vector of size seq_max with a clear configuration
    for(unsigned i = 0; i < length(row(align, 0)); ++i)  // maximum size of this string is length(mapline)
    {
        if (row(align, 0)[i]=='-') // In this choice the assumption that no alignment between - char can exist
        {
            ++gap0;
        }
        else if (row(align, 1)[i] == '-') // else if (row1[i]=='-')
        {
            ++gap1;
        }
        else
        {
            rnaAlign.mask[j] = std::make_pair(i + row0Begin - gap0, i + row1Begin - gap1);
            ++j;
        }
    }
    rnaAlign.maskIndex = j;
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

void computeUpperBound(TRnaAlign & rnaAlign)
{
    TScoreValue sum = 0;
    rnaAlign.slm = 0;
    for(unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        if (rnaAlign.upperBoundVect[i].maxProbScoreLine > 0)
        {
            sum += rnaAlign.upperBoundVect[i].maxProbScoreLine;
            if (rnaAlign.upperBoundVect[i].seq1Index !=
                rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq2IndexPairLine].seq1IndexPairLine)
            {
                ++rnaAlign.slm;
            }
        }
    }
    rnaAlign.upperBound = sum;
};

// ----------------------------------------------------------------------------
// Function computeBound()
// ----------------------------------------------------------------------------

void computeBound(TRnaAlign & rnaAlign)
{
    TScoreValue sumU = 0;
    TScoreValue sumL = 0;
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

void computeBounds(TRnaAlign & rnaAlign, TMapVect & lowerBound4Lemon)
{
    RnaStructureGraph & graph1 = rnaAlign.bppGraphH;
    RnaStructureGraph & graph2 = rnaAlign.bppGraphV;
    TScoreValue edgesProb, halfEdgesProb;
//  Clear the maxProbScoreLine of the upper bound
    for(unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i)
    {
        rnaAlign.upperBoundVect[i].maxProbScoreLine = 0;
    }
    for(unsigned i = 0; i < rnaAlign.maskIndex - 1; ++i)
    {
        String<unsigned> adjVect1;
        getVertexAdjacencyVector(adjVect1, graph1.inter, rnaAlign.mask[i].first);
        for (unsigned j = 0; j < length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)
        {
            for (unsigned w = i + 1; w < rnaAlign.maskIndex; ++w)
            {
                if (adjVect1[j] == rnaAlign.mask[w].first)
                {
                    String<unsigned> adjVect2;
                    getVertexAdjacencyVector(adjVect2, graph2.inter, rnaAlign.mask[i].second);
                    for (unsigned z = 0; z < length(adjVect2) && adjVect2[z] >= rnaAlign.mask[i].second; ++z)
                    {
                        if (adjVect2[z] == rnaAlign.mask[w].second)
                        {
                            edgesProb = (cargo(findEdge(graph1.inter, rnaAlign.mask[i].first, rnaAlign.mask[w].first)) +
                                         cargo(findEdge(graph2.inter, rnaAlign.mask[i].second, rnaAlign.mask[w].second)));
                            // Add edge for the LowerBoundGraph
                            lowerBound4Lemon[i][w] = edgesProb;
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
    computeUpperBound(rnaAlign);
}

// ----------------------------------------------------------------------------
// Function computeBounds() version that make use of an approximation of the MWM
// ----------------------------------------------------------------------------

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
    for (unsigned i = 0; i < rnaAlign.maskIndex - 1; ++i)
    {
        String<unsigned> adjVect1;
        getVertexAdjacencyVector(adjVect1, graph1.inter, rnaAlign.mask[i].first);
        for (unsigned j = 0; j < length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)
        {
            for (unsigned w = i + 1; w < rnaAlign.maskIndex; ++w)
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
    computeBound(rnaAlign);
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
    //TODO embed the LEMON GRAPH directly in this point in order to do not copy the structures if the lemon::MWM is chosen
    clearVertices(rnaAlign.lowerBoundGraph);
    // Add vertex for the LowerBoundGraph
    for(unsigned i = 0; i < rnaAlign.maskIndex; ++i)
    {
        addVertex(rnaAlign.lowerBoundGraph);
    }
    unsigned ll = 0;
    for (unsigned i = 0; i < rnaAlign.maskIndex - 1; ++i)
    {
        String<unsigned> adjVect1;
        getVertexAdjacencyVector(adjVect1, graph1.inter, rnaAlign.mask[i].first);
        for(unsigned j = 0; j< length(adjVect1) && adjVect1[j] >= rnaAlign.mask[i].first; ++j)
        {
            for (unsigned w = i + 1; w < rnaAlign.maskIndex; ++w)
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
                            std::cout << "Sum of probabilities = " << edgesProb << std::endl;

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

void saveBestAlignMinBound(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, unsigned index)
{
//    if ((rnaAlign.upperBound - rnaAlign.lowerBound) < (rnaAlign.upperMinBound - rnaAlign.lowerMinBound))
    if( rnaAlign.stepSize <= rnaAlign.stepSizeMinBound)  //TODO check if this <= is expensive
    {
        rnaAlign.itMinBounds = index; //to be used for the best lower bound
        rnaAlign.lowerMinBound = rnaAlign.lowerBound;
        rnaAlign.upperMinBound = rnaAlign.upperBound;
        rnaAlign.stepSizeMinBound = rnaAlign.stepSize;
        rnaAlign.bestAlignMinBounds = align;
        rnaAlign.bestAlignScoreMinBounds = alignScore;
    }
}

void updateLambda(TRnaAlign & rnaAlign) {
//    std::cout << "updateLambda function" << std::endl;
    for (unsigned i = 0; i < length(rnaAlign.upperBoundVect); ++i) {
        if (rnaAlign.upperBoundVect[i].maxProbScoreLine > 0) {
            // the edges are not paired
            if (rnaAlign.upperBoundVect[i].seq1Index !=
                rnaAlign.upperBoundVect[rnaAlign.upperBoundVect[i].seq1IndexPairLine].seq1IndexPairLine)
            {
                if (rnaAlign.upperBoundVect[i].seq1Index < rnaAlign.upperBoundVect[i].seq1IndexPairLine)
                {
// TODO check if this strategy is properly working a positive score is assigned to the left-side alignments.
// Maybe a double side strategy should be tested
                    rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].step += rnaAlign.stepSize;
                    // Note, the default initializer is callet the fist time that set the value to 0
                } else {
                    rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].step -= rnaAlign.stepSize;
                    // Note, the default initializer is callet the fist time that set the value to 0
                }
            }
// Save the maximum interaction weight to be used for the computation of profit of a line
            if (rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].maxProbScoreLine <
                rnaAlign.upperBoundVect[i].maxProbScoreLine)
            {
                rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].maxProbScoreLine
                        = rnaAlign.upperBoundVect[i].maxProbScoreLine;
                rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].seq1IndexPairLine
                        = rnaAlign.upperBoundVect[i].seq1IndexPairLine;
                rnaAlign.lamb[rnaAlign.upperBoundVect[i].seq1Index].map[i].seq2IndexPairLine
                        = rnaAlign.upperBoundVect[i].seq2IndexPairLine;
            }
        }
    }
}
#endif //_INCLUDE_ALIGNMENT_EDGES_H_
