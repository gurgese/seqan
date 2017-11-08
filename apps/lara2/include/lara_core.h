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
bool createMask(TRnaAlign & rnaAlign, TAlign const & align, TOptions const & options)
{
    typedef typename Iterator<Gaps<TSequence, seqan::ArrayGaps> const, Standard>::Type TGapsIter;

    Row<TAlign>::Type row0 = row(align, 0);
    Row<TAlign>::Type row1 = row(align, 1);

    // Get iterators.
    TGapsIter it0      = begin(row0);
    TGapsIter itEnd0   = end(row0);
    TGapsIter it1      = begin(row1);
    TGapsIter itEnd1   = end(row1);

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;
    // Keep track of the sequence positions for mask.
    std::pair<unsigned, unsigned> sourcePos(0u, 0u);
    // True, if this function changes the recent mask, False otherwise.
    bool changedMask = false;
    // Sum up the sequence and gap score.
    rnaAlign.sequenceScore = 0;

    for (unsigned lineCount = 0u; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1)
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
            ++sourcePos.second;
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
            ++sourcePos.first;
        }
        else
        {
            isGapOpen1 = false;
        }

        // Match or mismatch
        if (!isGap(it0) && !isGap(it1))
        {
            // create mask entry
            if (lineCount >= length(rnaAlign.mask))
            {
                appendValue(rnaAlign.mask, sourcePos);
                changedMask = true;
            }
            else if (rnaAlign.mask[lineCount] != sourcePos)
            {
                rnaAlign.mask[lineCount] = sourcePos;
                changedMask = true;
            }

            ++lineCount;
            ++sourcePos.first;
            ++sourcePos.second;
            rnaAlign.sequenceScore += score(options.laraScoreMatrix, *it0, *it1);
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);
    return changedMask;
}

// ----------------------------------------------------------------------------
// Function computeLowerBound()
// ----------------------------------------------------------------------------

void computeLowerBound(TRnaAlign & rnaAlign, TMapVect * lowerBound4Lemon)
{
    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;
    // iterate all lines that are present in the alignment
//    for (std::pair<unsigned, unsigned> const & line : rnaAlign.mask)
    //std::cout << graph1 << std::endl;
    //std::cout << graph2 << std::endl;
    String<unsigned> vectOut1, vectOut2;
//    unsigned nClosedLoops = 0;
//    unsigned nLinesInteraction = 0;
    for(unsigned i = 0; i < length(rnaAlign.mask) - 1; ++i)
    {
        seqan::String<std::pair<unsigned, unsigned> > seq2pos;
        //std::cout << rnaAlign.mask[i].first << " - " << rnaAlign.mask[i].second << std::endl;
        if (degree(graph1, rnaAlign.mask[i].first) > 0 && degree(graph2, rnaAlign.mask[i].second) > 0) {
            getVertexAdjacencyVector(vectOut1, graph1, rnaAlign.mask[i].first);
            for (unsigned x = 0; x < length(vectOut1); ++x) {
                for (unsigned j = i + 1; j < length(rnaAlign.mask); ++j) {
                    if (vectOut1[x] == rnaAlign.mask[j].first)
                    {
//                        std::cout << rnaAlign.mask[j].first << " : " << rnaAlign.mask[j].second << std::endl;
                        appendValue(seq2pos, std::make_pair(rnaAlign.mask[j].first, rnaAlign.mask[j].second));
                    }
                }
            }
//            ++nLinesInteraction;
        }
        if(length(seq2pos) > 0)
        {
            getVertexAdjacencyVector(vectOut2, graph2, rnaAlign.mask[i].second);
            for (unsigned w = 0; w < length(seq2pos); ++w)
            {
                for (unsigned y = 0; y < length(vectOut2); ++y)
                if (seq2pos[w].second == vectOut2[y])
                {
//                    std::cout << rnaAlign.mask[i].first << " - " << seq2pos[w].first << " = "
//                              << cargo(findEdge(graph1, rnaAlign.mask[i].first, seq2pos[w].first)) << " | "
//                              << rnaAlign.mask[i].second << " - " << seq2pos[w].second << " = "
//                              << cargo(findEdge(graph2, rnaAlign.mask[i].second, seq2pos[w].second))
//                              << std::endl;
//                    ++nClosedLoops;
                    (*lowerBound4Lemon)[rnaAlign.mask[i].first][seq2pos[w].first] = cargo(findEdge(graph1, rnaAlign.mask[i].first, seq2pos[w].first))
                                                                      + cargo(findEdge(graph2, rnaAlign.mask[i].second, seq2pos[w].second));

/*                    if (degree(graph1, rnaAlign.mask[j].first) > 0 && degree(graph2, rnaAlign.mask[j].second) > 0)
                    std::cout << "Interaction match: " << rnaAlign.mask[i].first+1 << " - " << rnaAlign.mask[i].second+1
                              << " | " << rnaAlign.mask[j].first+1 << " - " << rnaAlign.mask[j].second+1
                              << " // " << findEdge(graph1, rnaAlign.mask[i].first, rnaAlign.mask[j].first)
                              << " - " << findEdge(graph2, rnaAlign.mask[i].second, rnaAlign.mask[j].second)
                              << std::endl;
*/
                }
            }

        }
//        rnaAlign.slm = nLoops + 1 - (nClosedLoops * 2); // TODO check this value: It should be referred to the alignment line only or to all the lambda?
    }
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

//                std::cerr << "Interaction match: " << line.first+1 << " - " << line.second+1
//                          << " | " << value(adj_it1)+1 << " - " << value(adj_it2)+1 << "\tprob = "
//                          << edgeWeight1 << " + " << edgeWeight-edgeWeight1 << "\n";

                // for upper bound do not care if interactions are closed by a line
                if (rnaAlign.weightLineVect[line.second].maxProbScoreLine < edgeWeight / 2.0)
                {
                    rnaAlign.weightLineVect[line.second].maxProbScoreLine = edgeWeight / 2.0;
                    rnaAlign.weightLineVect[line.second].seq1Index = line.first;
                    rnaAlign.weightLineVect[line.second].seq1IndexPairLine = value(adj_it1);
                    rnaAlign.weightLineVect[line.second].seq2IndexPairLine = value(adj_it2);
//                    std::cerr << "updated\n";
                }
            }
        }
    }
}

void computeLowerBoundGreedy(TMapVect & interactions, TRnaAlign & rnaAlign) //TODO to be updated and verifyed
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
// Function saveBestAlignMinBound()
// ----------------------------------------------------------------------------

void saveBestAlignMinBound(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, unsigned index)
{
//    if ((rnaAlign.upperBound - rnaAlign.lowerBound) < (rnaAlign.upperMinBound - rnaAlign.lowerMinBound))
    if (rnaAlign.stepSize < rnaAlign.forMinBound.stepSizeBound)  //TODO check if this <= is expensive
    {
//        std::cerr << "update best min bound" << std::endl;
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

// ----------------------------------------------------------------------------
// Function saveBestAlignScore()
// ----------------------------------------------------------------------------

void saveBestAlignScore(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
//    if ((rnaAlign.upperBound - rnaAlign.lowerBound) < (rnaAlign.upperMinBound - rnaAlign.lowerMinBound))
    if (rnaAlign.forScore.bestAlignScore < alignScore)
    {
//        std::cerr << "update best score" << std::endl;
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

void saveBestAlignMinDiff(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
    if (rnaAlign.bestUpperBound - rnaAlign.bestLowerBound < rnaAlign.forMinDiff.upperBound - rnaAlign.forMinDiff.lowerBound)
    {
        rnaAlign.forMinDiff.it = index; //to be used for the best lower bound
        rnaAlign.forMinDiff.lowerBound = rnaAlign.bestLowerBound;
        rnaAlign.forMinDiff.upperBound = rnaAlign.bestUpperBound;
        rnaAlign.forMinDiff.stepSizeBound = rnaAlign.stepSize;
        rnaAlign.forMinDiff.bestAlign = align;
        rnaAlign.forMinDiff.bestAlignScore = alignScore;
        rnaAlign.forMinDiff.weightLineVect = rnaAlign.weightLineVect;
        rnaAlign.forMinDiff.mask = rnaAlign.mask;
    }
}

void saveBestAligns(TRnaAlign & rnaAlign, TAlign const & align, TScoreValue alignScore, int index)
{
    saveBestAlignMinDiff(rnaAlign, align, alignScore, index);
//    saveBestAlignMinBound(rnaAlign, align, alignScore, index);
//    saveBestAlignScore(rnaAlign, align, alignScore, index);
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

// ----------------------------------------------------------------------------
// Function createNewLambdaLines()
// ----------------------------------------------------------------------------

void createNewLambdaLines(TRnaAlign & rnaAlign, bool const & saveFoundInterPair)
{
    Graph<Undirected<double> > & graph1 = rnaAlign.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = rnaAlign.bppGraphV.inter;

    // iterate all lines that are present in the alignment
    for (std::pair<unsigned, unsigned> const & line : rnaAlign.mask)
    {
//        std::cerr << "graph degree " << degree(graph1, line.first) << " - " << degree(graph2, line.second) << std::endl;
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
            }
            if (proceed)
            {
                lambWeightStruct & lambda = rnaAlign.lamb[line.first].map[line.second];
                updateLambdaLine(lambda.maxProbScoreLine1, lambda.seq1IndexInter, graph1, line.first);
                lambda.seq1IndexPairLine = line.first;
//            std::cerr << lambda.seq1IndexPairLine << " : " << lambda.seq1IndexInter << " | " << lambda.maxProbScoreLine1
//                      << " - " << line.second << std::endl;
                updateLambdaLine(lambda.maxProbScoreLine2, lambda.seq2IndexInter, graph2, line.second);
                lambda.seq2IndexPairLine = line.second;
//            std::cerr << lambda.seq2IndexPairLine << " : " << lambda.seq2IndexInter << " | " << lambda.maxProbScoreLine2
//                      << " - " << line.first << std::endl;
                lambda.maxProbScoreLine = (lambda.maxProbScoreLine1 + lambda.maxProbScoreLine2) / 2;

//                std::cerr << lambda.seq1IndexPairLine << " : " << lambda.seq2IndexPairLine << " = "
//                          << lambda.maxProbScoreLine << " | " << lambda.seq1IndexInter << " : "
//                          << lambda.seq2IndexInter << " (" << lambda.fromUBPairing << ")" << std::endl;
                lambda.fromUBPairing = false;
                if (saveFoundInterPair)
                    // if the map is not empty means that the line have been already evaluated in a previous iteration
                {
                    bool proceed2 = false;
                    if(rnaAlign.lamb[lambda.seq1IndexInter].map.count(lambda.seq2IndexInter) == 0)
                    {
                        proceed2 = true;
                    } else if (rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].maxProbScoreLine < lambda.maxProbScoreLine)
                    {
                        proceed2 = true;
                    }
                    if (proceed2)
                    {
                        lambWeightStruct & lambdaPair = rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter];
                        lambdaPair.seq1IndexPairLine = lambda.seq1IndexInter;
                        lambdaPair.seq2IndexPairLine = lambda.seq2IndexInter;
                        lambdaPair.seq1IndexInter = lambda.seq1IndexPairLine;
                        lambdaPair.seq2IndexInter = lambda.seq2IndexPairLine;
                        lambdaPair.maxProbScoreLine = lambda.maxProbScoreLine;
                        lambdaPair.fromUBPairing = true;

//                        std::cerr << lambdaPair.seq1IndexPairLine << " : " << lambdaPair.seq2IndexPairLine << " = "
//                                  << lambdaPair.maxProbScoreLine << " | " << lambdaPair.seq1IndexInter << " : "
//                                  << lambdaPair.seq2IndexInter << " (" << lambdaPair.fromUBPairing << ")" << std::endl;
                    }
                }
            }
        }
    }
    std::cerr << "lambda num elements: ";
    for(unsigned i=0; i < length(rnaAlign.lamb); ++i)
        std::cerr << length(rnaAlign.lamb[i].map) << "\t";
    std::cerr << std::endl;

}

// ----------------------------------------------------------------------------
// Function computeSlm()
// ----------------------------------------------------------------------------

template <typename TMask>
void computeSlm(TRnaAlign & rnaAlign, TMask & listUnclosedLoopMask)
{
    rnaAlign.slm = 0;
    for (std::pair<unsigned, unsigned> const & lineL : rnaAlign.mask)
    {
        if (rnaAlign.lamb[lineL.first].map.count(lineL.second) > 0) {
            lambWeightStruct &lambdaL = rnaAlign.lamb[lineL.first].map[lineL.second];

            bool closedCircle = false;
            for (std::pair<unsigned, unsigned> const & lineM : rnaAlign.mask)
            {
                if (lambdaL.seq1IndexInter == lineM.first && lambdaL.seq2IndexInter == lineM.second)
                    closedCircle = true;
            }
            if (!closedCircle)
            {
                rnaAlign.slm += 2;
                appendValue(listUnclosedLoopMask, lineL);
            }
        }

/*
        if (rnaAlign.lamb[line.first].map.count(line.second) > 0)
        {
            lambWeightStruct & lambda = rnaAlign.lamb[line.first].map[line.second];
            if (lambda.seq1IndexPairLine !=
                rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].seq1IndexInter
                || lambda.seq2IndexPairLine !=
                   rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].seq2IndexInter)
            {
                ++rnaAlign.slm;

                std::cout << line.first << " : " << line.second << " = " << lambda.maxProbScoreLine << " | ";
                std::cout << lambda.seq1IndexInter << " : "
                          << lambda.seq2IndexInter << " = "
                          << rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].maxProbScoreLine
                          << std::endl;
                std::cout << rnaAlign.slm << std::endl;

                appendValue(listUnclosedLoopMask, line);
            }
        }
*/
    }
    std::cerr << "slm = " << rnaAlign.slm << std::endl;
}

// ----------------------------------------------------------------------------
// Function computeSlm()
// ----------------------------------------------------------------------------

template <typename TMask>
void updateLambdaStep(TRnaAlign & rnaAlign, TMask & listUnclosedLoopMask)
{
    for (std::pair<unsigned, unsigned> const &line : listUnclosedLoopMask)
    {
        if (rnaAlign.lamb[line.first].map.count(line.second) > 0)
        {
            lambWeightStruct & lambda = rnaAlign.lamb[line.first].map[line.second];
            if (true || lambda.seq1IndexPairLine !=
                rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].seq1IndexInter
                || lambda.seq2IndexPairLine !=
                   rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].seq2IndexInter)
            {
                lambda.step -= rnaAlign.stepSize;
                rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].step += rnaAlign.stepSize;
            }
/*
            std::cout << line.first << " : " << line.second << " = " << lambda.maxProbScoreLine << " -> "<< lambda.step << " | ";
            std::cout << lambda.seq1IndexInter << " : "
                      << lambda.seq2IndexInter << " = "
                      << rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].maxProbScoreLine << " -> "
                      << rnaAlign.lamb[lambda.seq1IndexInter].map[lambda.seq2IndexInter].step
                      << std::endl;
*/
        }
    }
}


//////// TRASH

// ----------------------------------------------------------------------------
// Function computeUpperBound()
// ----------------------------------------------------------------------------

/*
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
 */

void updateLambdaOld(TRnaAlign & rnaAlign)  // TODO REMOVE THIS FUNCTION
{
    for (size_t i = 0; i < length(rnaAlign.weightLineVect); ++i)
    {
        struct weightLineStruct const & ub = rnaAlign.weightLineVect[i]; //TODO revrite it with lambda only

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
