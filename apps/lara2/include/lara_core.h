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
// Function createInterLines()
// ----------------------------------------------------------------------------

void createInterLines(OutgoingInteractions & interactions,
                      NewLambdasVect & newLambdas,
                      int & numberOfEdgesInClosedLoops,
                      RnaInteractionGraph const & bppH,
                      RnaInteractionGraph const & bppV)
{
//    std::cout << bppH << std::endl;
//    std::cout << bppV << std::endl;
    for (unsigned ntSeq1Idx = 0; ntSeq1Idx < length(interactions); ++ntSeq1Idx)
    {
        //PositionPair const & lineL = traits.lines[lineLIdx];

        // skip if there are no interactions
        if (degree(bppH, ntSeq1Idx) == 0) // || degree(bppV, lineL.second) == 0)
            continue;
        // check for all adjacent pairs whether they are a line in the current alignment
        for (RnaAdjacencyIterator adjItH(bppH, ntSeq1Idx); !atEnd(adjItH); goNext(adjItH))
        {
            if(value(adjItH) > ntSeq1Idx) // only choose rightbound interactions
            {
                for (unsigned ntSeq2Idx = numVertices(bppV) - 1; ntSeq2Idx > 0; --ntSeq2Idx)
                {
                    // skip if there are no interactions
                    if (degree(bppV, ntSeq2Idx) == 0) // || degree(bppV, lineL.second) == 0)
                        continue;
                    for (RnaAdjacencyIterator adjItV(bppV, ntSeq2Idx); !atEnd(adjItV); goNext(adjItV))
                    {
                        if (value(adjItV) < ntSeq2Idx)
                        {
//                            std::cout  << ntSeq1Idx << " - ";
//                            std::cout << value(adjItH) << "\t";
//                            std::cout << "->" << ntSeq2Idx << "\t";
//                            std::cout << " - " << value(adjItV) << "\t || \t";
                            double strScore = ( cargo(findEdge(bppH, ntSeq1Idx, value(adjItH))) +
                                              cargo(findEdge(bppV, ntSeq2Idx, value(adjItV))) ) / 2;
                            InterLinePosWeight ilpw;

                            //Values to be used for th count of lines that do not close any loop
//                            if(interactions[ntSeq1Idx].count(value(adjItV)) == 0)
//                                ++numberOfEdgesInClosedLoops;
//                            if(interactions[value(adjItH)].count(ntSeq2Idx) == 0)
//                                ++numberOfEdgesInClosedLoops;

                            if(interactions[ntSeq1Idx].count(value(adjItV)) == 0 ||
                                    interactions[ntSeq1Idx][value(adjItV)].weight < strScore)
                            {
                                interactions[ntSeq1Idx][value(adjItV)].weight = strScore;
                                interactions[ntSeq1Idx][value(adjItV)].lineM = std::make_pair (value(adjItH), ntSeq2Idx);
                                //interactions[ntSeq1Idx][value(adjItV)].lineEnd.push_back(std::make_pair(value(adjItH), ntSeq2Idx));

                                newLambdas[ntSeq1Idx][value(adjItV)].posSeq1 = value(adjItH);
                                newLambdas[ntSeq1Idx][value(adjItV)].posSeq2 = ntSeq2Idx;
                                newLambdas[ntSeq1Idx][value(adjItV)].weight = strScore;
                            }
                            ilpw.posSeq1 = value(adjItH);
                            ilpw.posSeq2 = ntSeq2Idx;
                            ilpw.weight = strScore;
                            interactions[ntSeq1Idx][value(adjItV)].lineEnd[std::make_pair(value(adjItH), ntSeq2Idx)] = ilpw;

                            if(interactions[value(adjItH)].count(ntSeq2Idx) == 0 ||
                                    interactions[value(adjItH)][ntSeq2Idx].weight < strScore)
                            {
                                interactions[value(adjItH)][ntSeq2Idx].weight = strScore;
                                interactions[value(adjItH)][ntSeq2Idx].lineM = std::make_pair(ntSeq1Idx, value(adjItV));
                                //interactions[value(adjItH)][ntSeq2Idx].lineBegin.push_back(std::make_pair(ntSeq1Idx, value(adjItV)));

                                newLambdas[value(adjItH)][ntSeq2Idx].posSeq1 = ntSeq1Idx;
                                newLambdas[value(adjItH)][ntSeq2Idx].posSeq2 = value(adjItV);
                                newLambdas[value(adjItH)][ntSeq2Idx].weight = strScore;
                            }
                            ilpw.posSeq1 = ntSeq1Idx;
                            ilpw.posSeq2 = value(adjItV);
                            interactions[value(adjItH)][ntSeq2Idx].lineBegin[std::make_pair(ntSeq1Idx, value(adjItV))] = ilpw;
/*                            std::cout  << interactions[value(adjItH)][ntSeq2Idx].lineBegin.back().posSeq1 << " - ";
                            std::cout << interactions[value(adjItH)][ntSeq2Idx].lineBegin.back().posSeq2 << "\t";
                            std::cout << "->" << interactions[ntSeq1Idx][value(adjItV)].lineEnd.back().posSeq1 << "\t";
                            std::cout << " - " << interactions[ntSeq1Idx][value(adjItV)].lineEnd.back().posSeq2 << " = "
                                      << interactions[value(adjItH)][ntSeq2Idx].lineBegin.back().weight << "\t\t"
                                    << interactions[value(adjItH)][ntSeq2Idx].lineM.first << " - "
                                    << interactions[value(adjItH)][ntSeq2Idx].lineM.second << "\n";
*/
                        }
                    }
                }
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function evaluateLines()
// ----------------------------------------------------------------------------

bool evaluateLines(RnaAlignmentTraits & traits, RnaAlignment const & align, LaraOptions const & options)
{
    typedef typename Iterator<Gaps<Rna5String, ArrayGaps> const, Standard>::Type TGapsIter;

    Row<RnaAlignment>::Type row0 = row(align, 0);
    Row<RnaAlignment>::Type row1 = row(align, 1);

    // Get iterators.
    TGapsIter it0    = begin(row0);
    TGapsIter itEnd0 = end(row0);
    TGapsIter it1    = begin(row1);
    TGapsIter itEnd1 = end(row1);

    // State whether we have already opened a gap.
    bool isGapOpen0 = false;
    bool isGapOpen1 = false;

    // Keep track of the sequence positions for lines.
    PositionPair sourcePos(0u, 0u);

    // True, if this function changes the recent set of lines, False otherwise.
//    bool changedLines = false;
    // Sum up the sequence and gap score.
    traits.sequenceScore = 0.0;
    clear(traits.lines);

    for (unsigned lineCount = 0u; it0 != itEnd0 && it1 != itEnd1; ++it0, ++it1)
    {
        // Gaps in first sequence
        if (isGap(it0))
        {
            if (!isGapOpen0)
            {
                traits.sequenceScore += options.laraGapOpen;
            }
            else
            {
                traits.sequenceScore += options.laraGapExtend;
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
                traits.sequenceScore += options.laraGapOpen;
            }
            else
            {
                traits.sequenceScore += options.laraGapExtend;
            }
            isGapOpen1 = true;
            traits.pos2line[sourcePos.first] = UINT_MAX;
            ++sourcePos.first;
        }
        else
        {
            isGapOpen1 = false;
        }

        // Match or mismatch
        if (!isGap(it0) && !isGap(it1))
        {
            // create a line
            appendValue(traits.lines, sourcePos);
            traits.pos2line[sourcePos.first] = lineCount;
//            if (lineCount >= length(traits.lines))
//            {
//                appendValue(traits.lines, sourcePos);
//                changedLines = true;
//            }
//            else if (traits.lines[lineCount] != sourcePos)
//            {
//                traits.lines[lineCount] = sourcePos;
//                changedLines = true;
//            }
            traits.sequenceScore += score(traits.structureScore.matrix, *it0, *it1);

            ++lineCount;
            ++sourcePos.first;
            ++sourcePos.second;
//            traits.sequenceScore += score(options.laraScoreMatrix, *it0, *it1);
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);
//    return changedLines;
    return true;
}

// ----------------------------------------------------------------------------
// Function computeLowerBound()
// ----------------------------------------------------------------------------

double computeLowerBoundGreedy(RnaAlignmentTraits & traits)
{
    // add vertices to interaction match graph
    clear(traits.interactionMatchingGraph);
    forEach(traits.lines, [&traits] (PositionPair const &) { addVertex(traits.interactionMatchingGraph); });

    RnaInteractionGraph const & bppH = traits.bppGraphH.inter;
    RnaInteractionGraph const & bppV = traits.bppGraphV.inter;

    traits.numberOfSubgradients = 0;

    //TODO modify this piece of code for keeping the already computed closed loops from traits.interactions
    // add edges: go through all lines
    for (unsigned lineLIdx = 0; lineLIdx < length(traits.lines); ++lineLIdx)
    {
        PositionPair const & lineL = traits.lines[lineLIdx];

        // skip if there are no interactions
        if (degree(bppH, lineL.first) == 0 || degree(bppV, lineL.second) == 0)
            continue;

        // count number of outgoing max interctions
        traits.numberOfSubgradients += 2;

        // check for all adjacent pairs whether they are a line in the current alignment
        for (RnaAdjacencyIterator adjItH(bppH, lineL.first); !atEnd(adjItH); goNext(adjItH))
        {
            // search for an existing lineM that builds a closed cycle with lineL
            unsigned lineMIdx = traits.pos2line[value(adjItH)];
            if (lineMIdx != UINT_MAX && lineLIdx < lineMIdx)
            {
                PositionPair const & lineM = traits.lines[lineMIdx];
                for (RnaAdjacencyIterator adjItV(bppV, lineL.second); !atEnd(adjItV); goNext(adjItV))
                {
                    if (lineM.second == value(adjItV)) // found a closed cycle
                    {
                        addEdge(traits.interactionMatchingGraph, lineLIdx, lineMIdx,
                                cargo(findEdge(bppH, lineL.first, lineM.first)) +
                                cargo(findEdge(bppV, lineL.second, lineM.second)));
                        // break out of the two inner for loops and continue with next lineL
                        //goto nextLineL;
                    }
                }
            }
        }

//        nextLineL:
//        continue;
    }

    double mwmScore = maximumWeightedMatchingGreedy<5>(traits.isInMwmSolution, traits.interactionMatchingGraph);
    for (bool edge : traits.isInMwmSolution)
        if(edge)
            traits.numberOfSubgradients -= 2;
    SEQAN_ASSERT_GEQ(traits.numberOfSubgradients, 0);

    return mwmScore;
}

// ----------------------------------------------------------------------------
// Function saveBestAlignMinBound()
// ----------------------------------------------------------------------------

void saveBestAlignMinBound(RnaAlignmentTraits & traits, RnaAlignment const & align, double alignScore, unsigned index)
{
//    if ((traits.upperBound - traits.lowerBound) < (traits.upperMinBound - traits.lowerMinBound))
    if (traits.stepSize < traits.forMinBound.stepSizeBound)  //TODO check if this <= is expensive
    {
//        std::cerr << "update best min bound" << std::endl;
        traits.forMinBound.it = index; //to be used for the best lower bound
        traits.forMinBound.lowerBound = traits.lowerBound;
        traits.forMinBound.upperBound = traits.upperBound;
        traits.forMinBound.stepSizeBound = traits.stepSize;
        traits.forMinBound.bestAlign = align;
        traits.forMinBound.bestAlignScore = alignScore;
//        traits.forMinBound.weightLineVect = traits.weightLineVect;
        traits.forMinBound.lines = traits.lines;
    }
}

// ----------------------------------------------------------------------------
// Function saveBestAlignScore()
// ----------------------------------------------------------------------------
/*
void saveBestAlignScore(RnaAlignmentTraits & traits, RnaAlignment const & align, double alignScore, int index)
{
//    if ((traits.upperBound - traits.lowerBound) < (traits.upperMinBound - traits.lowerMinBound))
    if (traits.forScore.bestAlignScore < alignScore)
    {
//        std::cerr << "update best score" << std::endl;
        traits.forScore.it = index; //to be used for the best lower bound
        traits.forScore.lowerBound = traits.lowerBound;
        traits.forScore.upperBound = traits.upperBound;
        traits.forScore.stepSizeBound = traits.stepSize;
        traits.forScore.bestAlign = align;
        traits.forScore.bestAlignScore = alignScore;
//        traits.forScore.weightLineVect = traits.weightLineVect;
        traits.forScore.lines = traits.lines;
    }
}
 */

void saveBestAlignMinDiff(RnaAlignmentTraits & traits, RnaAlignment const & align, double alignScore, int index)
{
    if (traits.bestUpperBound - traits.bestLowerBound < traits.forMinDiff.upperBound - traits.forMinDiff.lowerBound)
    {
        traits.forMinDiff.it = index; //to be used for the best lower bound
        traits.forMinDiff.lowerBound = traits.bestLowerBound;
        traits.forMinDiff.upperBound = traits.bestUpperBound;
        traits.forMinDiff.stepSizeBound = traits.stepSize;
        traits.forMinDiff.bestAlign = align;
        traits.forMinDiff.bestAlignScore = alignScore;
//        traits.forMinDiff.weightLineVect = traits.weightLineVect;
        traits.forMinDiff.lines = traits.lines;
    }
}

void saveBestAligns(RnaAlignmentTraits & traits, RnaAlignment const & align, double alignScore, int index)
{
    saveBestAlignMinDiff(traits, align, alignScore, index);
    saveBestAlignMinBound(traits, align, alignScore, index);
//    saveBestAlignScore(traits, align, alignScore, index);
}

void findMaxWeight(double & maxWeight, unsigned & partnerIndex, RnaInteractionGraph const & graph, unsigned position)
{
    for (RnaAdjacencyIterator adj_it(graph, position); !atEnd(adj_it); goNext(adj_it))
    {
        double edgeWeight = cargo(findEdge(graph, position, value(adj_it)));
        if (maxWeight < edgeWeight)
        {
            maxWeight = edgeWeight;
            partnerIndex = value(adj_it);
        }
    }
}


// ----------------------------------------------------------------------------
// Function evaluateInteractions()
// ----------------------------------------------------------------------------

void evaluateInteractions(RnaAlignmentTraits & traits, unsigned const & iter)
{
    Graph<Undirected<double> > & graph1 = traits.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = traits.bppGraphV.inter;

    //bool flagUpdateClosedLoops = false;
    // iterate all lines that are present in the alignment
    for (PositionPair const & line : traits.lines)
    {
//        std::cerr << "graph degree " << degree(graph1, line.first) << " - " << degree(graph2, line.second) << std::endl;
        // outgoing interactions in the first sequence

        if (degree(graph1, line.first) == 0 || degree(graph2, line.second) == 0)
            continue;

        if (traits.interactions[line.first].count(line.second) > 0          // skip if interaction exists already
            && !traits.interactions[line.first][line.second].fromUBPairing) // and was evaluated before
            continue;

        InterLine & interaction = traits.interactions[line.first][line.second];
        interaction.lineL = line;
        interaction.iterUpdate = iter;
        findMaxWeight(interaction.maxProbScoreLine1, interaction.lineM.first, graph1, line.first);

//            std::cerr << interaction.lineL.first << " : " << interaction.lineM.first << " | " << interaction.maxProbScoreLine1
//                      << " - " << line.second << std::endl;
        findMaxWeight(interaction.maxProbScoreLine2, interaction.lineM.second, graph2, line.second);

//            std::cerr << interaction.lineL.second << " : " << interaction.lineM.second << " | " << interaction.maxProbScoreLine2
//                      << " - " << line.first << std::endl;
        interaction.weight = (interaction.maxProbScoreLine1 + interaction.maxProbScoreLine2) / 2;

//                std::cerr << interaction.lineL.first << " : " << interaction.lineL.second << " = "
//                          << interaction.weight << " | " << interaction.lineM.first << " : "
//                          << interaction.lineM.second << " (" << interaction.fromUBPairing << ")" << std::endl;
        interaction.fromUBPairing = false;
        interaction.lambdaValue = 0.0;
        //flagUpdateClosedLoops = true;


        if (traits.interactions[interaction.lineM.first].count(interaction.lineM.second) == 0
            || traits.interactions[interaction.lineM.first][interaction.lineM.second].weight < interaction.weight)
        {
            InterLine & pair = traits.interactions[interaction.lineM.first][interaction.lineM.second];
            pair.lineL = interaction.lineM;
            pair.lineM = interaction.lineL;
            pair.weight = interaction.weight;
            pair.fromUBPairing = true;
            pair.lambdaValue = 0.0;
            pair.iterUpdate = iter;
        }
    }
/*
    // update which loops are closed
    for (PositionPair const & line : traits.lines)
    {
        if (traits.interactions[line.first].count(line.second) == 0)
            continue;

        RnaInteraction & interaction = traits.interactions[line.first][line.second];
        SEQAN_ASSERT(traits.interactions[interaction.lineM.first].count(interaction.lineM.second) > 0);
        RnaInteraction & pairedInteraction = traits.interactions[interaction.lineM.first][interaction.lineM.second];
        if (interaction.lineL == pairedInteraction.lineM)
        {
            interaction.closedLoop = true;
            pairedInteraction.closedLoop = true;
        }
        else
        {
            interaction.closedLoop = false;
            pairedInteraction.closedLoop = false;
        }

    }
    */
    if (iter == 0)
        return;

    std::cerr << "upper bound contributions:" << std::endl;
    double sum = 0;
    double sumLam = 0;
    double sumSol = 0;

    for (unsigned lineLIdx = 0u; lineLIdx < length(traits.lines); ++lineLIdx)
    {
        PositionPair const & line = traits.lines[lineLIdx];
        if (traits.interactions[line.first].count(line.second) == 0)
        {
            ++lineLIdx;
            continue;
        }

        InterLine & interaction = traits.interactions[line.first][line.second];
        sum += interaction.weight;
        sumLam += interaction.lambdaValue;

        bool inSolution = false;
        unsigned lineMIdx = traits.pos2line[interaction.lineM.first];
        if (lineMIdx != UINT_MAX && interaction.lineM == traits.lines[lineMIdx])
        {
            // line M exists in current alignment
            PositionPair const & lineM = traits.lines[lineMIdx];
            if (traits.interactions[lineM.first][lineM.second].lineM == line) // closed cycle of maximum weights
            {
                auto const & edgeLM = findEdge(traits.interactionMatchingGraph, lineLIdx, lineMIdx);
                if (edgeLM != 0 && traits.isInMwmSolution[edgeLM->data_id]) // cycle is represented in valid solution
                    inSolution = true;
            }
        }

        if (inSolution)
            sumSol += interaction.lambdaValue;

        std::cerr << "max_interaction weight = " << interaction.weight
                  << " lambda = " << interaction.lambdaValue << " ("
                  << interaction.lineL.first << "," << interaction.lineL.second << " - "
                  << interaction.lineM.first << "," << interaction.lineM.second << ") "
                  << "\tin solution: " << (inSolution ? "yes" : "no") << std::endl;
    }

    std::cerr << "sum upper bound cargo  = " << sum << std::endl;
    std::cerr << "sum upper bound lambda = " << sumLam << std::endl;
    std::cerr << "sum lambda in solution = " << sumSol << std::endl;
}

// ----------------------------------------------------------------------------
// Function computeNumberOfSubgradients()
// ----------------------------------------------------------------------------

//void computeNumberOfSubgradients2(RnaAlignmentTraits & traits)
//{
//    traits.numberOfSubgradients = 0; //TODO verify what happens if no structure is provided this value to 0 can give problems in the division
//    for (unsigned idx = 0u; idx < length(traits.interactions); ++idx)
//    {
//        for (const auto &interPair : traits.interactions[idx])
//        {
//            std::cout << idx << " - " << interPair.first << " || ";
//            std::cout << interPair.second.lineM.first << " - "
//                      << interPair.second.lineM.second << "\t:\t";
//            std::cout << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineM.first << " - "
//                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineM.second << "\n";
//            traits.numberOfSubgradients += length(interPair.second.lineBegin) + length(interPair.second.lineEnd);
//            if (idx == traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineM.first &&
//                interPair.first == traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineM.second)
//            {
//                --traits.numberOfSubgradients;
//            }
//        }
//    }
//    std::cerr << "number of subgradients = " << traits.numberOfSubgradients << std::endl;
//}


// ----------------------------------------------------------------------------
// Function computeNumberOfSubgradients()
// ----------------------------------------------------------------------------

void computeNumberOfSubgradients(RnaAlignmentTraits & traits)
{
    traits.numberOfSubgradients = 0;
    for (unsigned lineLIdx = 0u; lineLIdx < length(traits.lines); ++lineLIdx) {
        PositionPair const & line = traits.lines[lineLIdx];
        if (traits.interactions[line.first].count(line.second) == 0)
            continue;

        InterLine & interaction = traits.interactions[line.first][line.second];
        unsigned lineMIdx = traits.pos2line[interaction.lineM.first];

        if (lineMIdx == UINT_MAX || interaction.lineM != traits.lines[lineMIdx])
        {
            // lineM is not contained in current alignment
            traits.numberOfSubgradients += 2; // also consider reverse direction (M -> L)
            //appendValue(unclosedLoops, lineLIdx);
            //std::cerr << "unclosed: (" << traits.lines[lineLIdx].first << "," << traits.lines[lineLIdx].second << ")\n";
        }
        else
        {
            auto const & edgeLM = findEdge(traits.interactionMatchingGraph, lineLIdx, lineMIdx);
            if (edgeLM == 0 || !traits.isInMwmSolution[edgeLM->data_id])
            {
                // lineM exists, but is not contained in the MWM solution
                traits.numberOfSubgradients += 1; // reverse direction will be processed in separate iteration
                //appendValue(unclosedLoops, lineLIdx);
            }
        }
    }
    std::cerr << "number of subgradients = " << traits.numberOfSubgradients << std::endl;
}
/*
        bool stucturalLine = false;

        if (interaction.closedLoop)
        {
            if (interaction.iterUpdate >= 0 &&
                interaction.iterUpdate == traits.interactions[interaction.lineM.first][interaction.lineM.second].iterUpdate)
            {
                if (interaction.lineL.first !=
                    traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.first
                    || interaction.lineL.second !=
                       traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.second)
                {
                    stucturalLine = true;
                }
            }
            else
            {
                stucturalLine = true;
            }
        }
        else
        {
            stucturalLine = true;
        }
        if (stucturalLine)
        {
            traits.numberOfSubgradients += 2; // 1; TODO Verify this value
            appendValue(unclosedLoops, line);
        }

    }
*/


// ----------------------------------------------------------------------------
// Function updateLambdaValues()
// ----------------------------------------------------------------------------

void updateLambdaValues(RnaAlignmentTraits & traits)
{
//    for (unsigned lineLIdx : unclosedLoops)
//    {
//        PositionPair const & lineL = traits.lines[lineLIdx];
//        if (traits.interactions[lineL.first].count(lineL.second) == 0)
//            continue;
//
//        RnaInteraction & interaction = traits.interactions[lineL.first][lineL.second];
//        interaction.lambdaValue -= traits.stepSize;
//
//        // Check whether lambda of paired lineM directs back to lineL
//        RnaInteraction & paired = traits.interactions[interaction.lineM.first][interaction.lineM.second];
//        if (paired.lineM == interaction.lineL)
//        {
//            paired.lambdaValue += traits.stepSize;
//        }
//    }

    for (unsigned lineLIdx = 0u; lineLIdx < length(traits.lines); ++lineLIdx)
    {
        PositionPair const & lineL = traits.lines[lineLIdx];
        if (traits.interactions[lineL.first].count(lineL.second) == 0u)
        {
            SEQAN_ASSERT_EQ(degree(traits.interactionMatchingGraph, lineLIdx), 0u);
            continue;
        }

        InterLine & interaction = traits.interactions[lineL.first][lineL.second];
        double lambdaChange = std::min(interaction.weight + interaction.lambdaValue, traits.stepSize);
        unsigned lineMIdx = traits.pos2line[interaction.lineM.first];

        if (lineMIdx != UINT_MAX && interaction.lineM == traits.lines[lineMIdx])
        {
            // line M exists in current alignment
            PositionPair const & lineM = traits.lines[lineMIdx];
            if (traits.interactions[lineM.first][lineM.second].lineM == lineL) // closed cycle of maximum weights
            {
                auto const & edgeLM = findEdge(traits.interactionMatchingGraph, lineLIdx, lineMIdx);
                SEQAN_ASSERT(edgeLM != 0);
                if (traits.isInMwmSolution[edgeLM->data_id]) // cycle is represented in valid solution
                {
                    interaction.lambdaValue = 0.0;
                    traits.interactions[interaction.lineM.first][interaction.lineM.second].lambdaValue = 0.0;
                    continue;
                }
            }

            for (RnaAdjacencyIterator adjIt(traits.interactionMatchingGraph, lineLIdx); !atEnd(adjIt); goNext(adjIt))
            {
                auto const & edgeLM = findEdge(traits.interactionMatchingGraph, lineLIdx, value(adjIt));
                if (traits.isInMwmSolution[edgeLM->data_id])
                {
                    SEQAN_ASSERT_LEQ(getCargo(edgeLM) / 2.0, interaction.weight);
                    lambdaChange = std::min(lambdaChange, interaction.weight + interaction.lambdaValue - (getCargo(edgeLM) / 2.0));
                }
            }
        }

        interaction.lambdaValue -= lambdaChange;
        traits.interactions[interaction.lineM.first][interaction.lineM.second].lambdaValue += lambdaChange;
    }
}

// ----------------------------------------------------------------------------
// Function updateLambdaValues()
// ----------------------------------------------------------------------------

void updateLambdaValues2(RnaAlignmentTraits & traits, LaraOptions const & options)
{
//    RnaInteractionGraph const & bppH = traits.bppGraphH.inter;
//    RnaInteractionGraph const & bppV = traits.bppGraphV.inter;
//    std::cout << "Evaluate Lambda " << std::endl;
    //std::vector<InterLinePosWeight> newLambdas;
    //InterLinePosWeight newLambdaElem; // In this structure the weight field is used for storing the best Lambda.
    //for (unsigned ntSeq1Idx = 0u; ntSeq1Idx < length(traits.interactions); ++ntSeq1Idx) // position in seq1
    for (std::map<unsigned, InterLine> & interaction : traits.interactions)
    {
        for (std::pair<unsigned const, InterLine> & interPair : interaction) // the current line [seq2 pos, InterLine]
        {
            InterLinePosWeight newLambdaElem; // In this structure the weight field is used for storing the best Lambda.
            for (std::pair<PositionPair const, InterLinePosWeight> & endPair : interPair.second.lineEnd) // scan pairs on the right
            {
                /*std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                          <<  "\t L2:" << endPair.second.posSeq1 << "-" << endPair.second.posSeq2 << " = WL:"
                          << endPair.second.weight << " LA2:"
                          << traits.interactions[endPair.second.posSeq1][endPair.second.posSeq2].lambdaValue
                          << std::endl;
                */
                // if uninitialised
                if( newLambdaElem.posSeq1 < 0 || ( (newLambdaElem.weight + newLambdaElem.lambdaValue <
                                                    endPair.second.weight + endPair.second.lambdaValue) ) )
                {
                    newLambdaElem.posSeq1 = endPair.second.posSeq1;
                    newLambdaElem.posSeq2 = endPair.second.posSeq2;
                    newLambdaElem.weight = endPair.second.weight;
                    newLambdaElem.lambdaValue = endPair.second.lambdaValue;
                }
            }
            for (std::pair<PositionPair const, InterLinePosWeight> & beginPair : interPair.second.lineBegin) // scan pairs on the left
            {
                /*std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                          <<  "\t L2:" << beginPair.second.posSeq1 << "-" << beginPair.second.posSeq2 << " = WL:"
                          << beginPair.second.weight << " LA2:"
                          << traits.interactions[beginPair.second.posSeq1][beginPair.second.posSeq2].lambdaValue
                          << std::endl;
                */
                if( newLambdaElem.posSeq1 < 0 || ( (newLambdaElem.weight + newLambdaElem.lambdaValue <
                                                    beginPair.second.weight + beginPair.second.lambdaValue) ) )
                {
                    newLambdaElem.posSeq1 = beginPair.second.posSeq1;
                    newLambdaElem.posSeq2 = beginPair.second.posSeq2;
                    newLambdaElem.weight = beginPair.second.weight;
                    newLambdaElem.lambdaValue = beginPair.second.lambdaValue;
                }
            }
            interPair.second.weight = newLambdaElem.weight;
            interPair.second.lineM = std::make_pair (newLambdaElem.posSeq1, newLambdaElem.posSeq2);
            //interPair.second.lambdaValue = newLambdaElem.lambdaValue;
        }
        /*
        for (auto &interPair : traits.interactions[ntSeq1Idx])
        {
            std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                          <<  "\t L2:" << interPair.second.lineM.first << "-" << interPair.second.lineM.second << " = WL:"
                          << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].weight << " LA2:"
                          << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue
                         << std::endl;
            interPair.second.lambdaValue -= traits.stepSize;
            traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue += traits.stepSize;
            if(ntSeq1Idx < interPair.second.lineM.first)
            {
                interPair.second.lineEnd[std::make_pair (interPair.second.lineM.first, interPair.second.lineM.second)].lambdaValue += traits.stepSize;
                traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineBegin[std::make_pair (ntSeq1Idx, interPair.first)].lambdaValue -= traits.stepSize;
            } else
            {
                interPair.second.lineBegin[std::make_pair (interPair.second.lineM.first, interPair.second.lineM.second)].lambdaValue += traits.stepSize;
                traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineEnd[std::make_pair (ntSeq1Idx, interPair.first)].lambdaValue -= traits.stepSize;
            }
            std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                      <<  "\t L2:" << interPair.second.lineM.first << "-" << interPair.second.lineM.second << " = WL:"
                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].weight << " LA2:"
                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue
                      << std::endl;
        }
         */
    }
//    std::cout << "Update Lambda " << std::endl;
    for (unsigned ntSeq1Idx = 0u; ntSeq1Idx < length(traits.interactions); ++ntSeq1Idx)
    {
        for (std::pair<unsigned const, InterLine> & interPair : traits.interactions[ntSeq1Idx])
        {
//            std::cout << "before) L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue
//                      << " W1:" << interPair.second.weight
//                      << "\t L2:" << interPair.second.lineM.first << "-" << interPair.second.lineM.second << " LA2:"
//                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue
//                      << " = W2:"
//                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].weight
//                      << std::endl;

            interPair.second.lambdaValue -= traits.stepSize; // gamma
            traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue += traits.stepSize;
            if(ntSeq1Idx < interPair.second.lineM.first)
            {
                interPair.second.lineEnd[std::make_pair (interPair.second.lineM.first, interPair.second.lineM.second)].lambdaValue += traits.stepSize;
                traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineBegin[std::make_pair ((int)ntSeq1Idx, interPair.first)].lambdaValue -= traits.stepSize;
            } else
            {
                SEQAN_ASSERT(ntSeq1Idx != interPair.second.lineM.first);
                interPair.second.lineBegin[std::make_pair (interPair.second.lineM.first, interPair.second.lineM.second)].lambdaValue += traits.stepSize;
                traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lineEnd[std::make_pair ((int)ntSeq1Idx, interPair.first)].lambdaValue -= traits.stepSize;
            }
//            std::cout << "after)  L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
//                      <<  "\t L2:" << interPair.second.lineM.first << "-" << interPair.second.lineM.second << " LA2:"
//                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue
//                      << " = W2:"
//                      << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].weight
//                      << std::endl;

        }
        /*
        for (auto &interPair : traits.interactions[ntSeq1Idx])
        {
            InterLinePosWeight newLambdaElem; // In this structure the weight field is used for storing the best Lambda.
            for (auto &endPair : interPair.second.lineEnd)
            {
                std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                          <<  "\t L2:" << endPair.second.posSeq1 << "-" << endPair.second.posSeq2 << " = WL:"
                          << endPair.second.weight << " LA2:" //<< " = " << strScore << "\t";
                          << traits.interactions[endPair.second.posSeq1][endPair.second.posSeq2].lambdaValue
                          //<< " = " << strScore
                          << std::endl;
                totLamb += endPair.second.lambdaValue;
                if( newLambdaElem.posSeq1 < 0 || ( (newLambdaElem.weight + newLambdaElem.lambdaValue <
                        endPair.second.weight + endPair.second.lambdaValue) ) )
                {
                    newLambdaElem.posSeq1 = endPair.second.posSeq1;
                    newLambdaElem.posSeq2 = endPair.second.posSeq2;
                    newLambdaElem.weight = endPair.second.weight;
                    newLambdaElem.lambdaValue = endPair.second.lambdaValue;
                }
            }
            for (auto &beginPair : interPair.second.lineBegin)
            {
                std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue  << " W1:" << interPair.second.weight
                          <<  "\t L2:" << beginPair.second.posSeq1 << "-" << beginPair.second.posSeq2 << " = WL:"
                          << beginPair.second.weight << " LA2:" //<< " = " << strScore << "\t";
                          << traits.interactions[beginPair.second.posSeq1][beginPair.second.posSeq2].lambdaValue
                          //<< " = " << strScore
                          << std::endl;
                totLamb += beginPair.second.lambdaValue;
                if( newLambdaElem.posSeq1 < 0 || ( (newLambdaElem.weight + newLambdaElem.lambdaValue <
                        beginPair.second.weight + beginPair.second.lambdaValue) ) )
                {
                    newLambdaElem.posSeq1 = beginPair.second.posSeq1;
                    newLambdaElem.posSeq2 = beginPair.second.posSeq2;
                    newLambdaElem.weight = beginPair.second.weight;
                    newLambdaElem.lambdaValue = beginPair.second.lambdaValue;
                }
            }
            interPair.second.weight = newLambdaElem.weight;
            interPair.second.lineM = std::make_pair (newLambdaElem.posSeq1, newLambdaElem.posSeq2);
            interPair.second.lambdaValue = newLambdaElem.lambdaValue;
            usedLamb += interPair.second.lambdaValue;
        }
         */
    }
    if(options.verbose > 2)
    {
        double totLamb = 0;
//        std::cout << "Print new Lambda " << std::endl;
        for (unsigned ntSeq1Idx = 0u; ntSeq1Idx < length(traits.interactions); ++ntSeq1Idx)
        {
            for (auto &interPair : traits.interactions[ntSeq1Idx])
            {
//                for (auto &endPair : interPair.second.lineEnd)
//                {
//                    std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue
//                              << " W1:" << interPair.second.weight
//                              << "\t L2:" << endPair.second.posSeq1 << "-" << endPair.second.posSeq2 << " = WL:"
//                              << endPair.second.weight << " LA2:"
//                              << traits.interactions[endPair.second.posSeq1][endPair.second.posSeq2].lambdaValue
//                              << std::endl;
//                }
//                for (auto &beginPair : interPair.second.lineBegin)
//                {
//                    std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue
//                              << " W1:" << interPair.second.weight
//                              << "\t L2:" << beginPair.second.posSeq1 << "-" << beginPair.second.posSeq2 << " = WL:"
//                              << beginPair.second.weight << " LA2:"
//                              << traits.interactions[beginPair.second.posSeq1][beginPair.second.posSeq2].lambdaValue
//                              << std::endl;
//                }
//                std::cout << "L1:" << ntSeq1Idx << "-" << interPair.first << " LA1:" << interPair.second.lambdaValue
//                              << " W1:" << interPair.second.weight
//                              << "\t L2:" << interPair.second.lineM.first << "-" << interPair.second.lineM.second << " = LA2:"
//                              << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].lambdaValue
//                              << " W2:"
//                              << traits.interactions[interPair.second.lineM.first][interPair.second.lineM.second].weight
//                              << std::endl;
                totLamb = totLamb + interPair.second.lambdaValue;
            }
        }
//        std::cout <<"Stepsize = " << traits.stepSize << " total lambda = " << totLamb << std::endl;
    }


}

#endif //_INCLUDE_ALIGNMENT_EDGES_H_
