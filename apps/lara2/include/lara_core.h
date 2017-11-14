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
// Function evaluateLines()
// ----------------------------------------------------------------------------

bool evaluateLines(RnaAlignmentTraits & traits, RnaAlignment const & align, LaraOptions const & options)
{
    typedef typename Iterator<Gaps<Rna5String, ArrayGaps> const, Standard>::Type TGapsIter;

    Row<RnaAlignment>::Type row0 = row(align, 0);
    Row<RnaAlignment>::Type row1 = row(align, 1);

    // Get iterators.
    TGapsIter it0      = begin(row0);
    TGapsIter itEnd0   = end(row0);
    TGapsIter it1      = begin(row1);
    TGapsIter itEnd1   = end(row1);

    // State whether we have already opened a gap.
    bool isGapOpen0 = false, isGapOpen1 = false;
    // Keep track of the sequence positions for lines.
    PositionPair sourcePos(0u, 0u);
    // True, if this function changes the recent set of lines, False otherwise.
    bool changedLines = false;
    // Sum up the sequence and gap score.
    traits.sequenceScore = 0;

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
            if (lineCount >= length(traits.lines))
            {
                appendValue(traits.lines, sourcePos);
                changedLines = true;
            }
            else if (traits.lines[lineCount] != sourcePos)
            {
                traits.lines[lineCount] = sourcePos;
                changedLines = true;
            }

            ++lineCount;
            ++sourcePos.first;
            ++sourcePos.second;
            traits.sequenceScore += score(options.laraScoreMatrix, *it0, *it1);
        }
    }
    SEQAN_ASSERT(it0 == itEnd0);
    SEQAN_ASSERT(it1 == itEnd1);
    return changedLines;
}

// ----------------------------------------------------------------------------
// Function computeLowerBound()
// ----------------------------------------------------------------------------

void prepareLowerBoundScores(InteractionScoreMap & validInteractions, RnaAlignmentTraits const & traits)
{
    RnaInteractionGraph const & graph1 = traits.bppGraphH.inter;
    RnaInteractionGraph const & graph2 = traits.bppGraphV.inter;
    // iterate all lines that are present in the alignment
//    for (PositionPair const & line : traits.lines)
    //std::cout << graph1 << std::endl;
    //std::cout << graph2 << std::endl;
    String<unsigned> vectOut1, vectOut2;
//    unsigned nClosedLoops = 0;
//    unsigned nLinesInteraction = 0;
    for(unsigned idx = 0; idx < length(traits.lines) - 1; ++idx)
    {
        seqan::String<PositionPair > seq2pos;
        //std::cout << traits.lines[idx].first << " - " << traits.lines[idx].second << std::endl;
        if (degree(graph1, traits.lines[idx].first) > 0 && degree(graph2, traits.lines[idx].second) > 0) {
            getVertexAdjacencyVector(vectOut1, graph1, traits.lines[idx].first);
            for (unsigned x = 0; x < length(vectOut1); ++x) {
                for (unsigned j = idx + 1; j < length(traits.lines); ++j) {
                    if (vectOut1[x] == traits.lines[j].first)
                    {
//                        std::cout << traits.lines[j].first << " : " << traits.lines[j].second << std::endl;
                        appendValue(seq2pos, std::make_pair(traits.lines[j].first, traits.lines[j].second));
                    }
                }
            }
//            ++nLinesInteraction;
        }
        if(length(seq2pos) > 0)
        {
            getVertexAdjacencyVector(vectOut2, graph2, traits.lines[idx].second);
            for (unsigned w = 0; w < length(seq2pos); ++w)
            {
                for (unsigned y = 0; y < length(vectOut2); ++y)
                if (seq2pos[w].second == vectOut2[y])
                {
//                    std::cout << traits.lines[idx].first << " - " << seq2pos[w].first << " = "
//                              << cargo(findEdge(graph1, traits.lines[idx].first, seq2pos[w].first)) << " | "
//                              << traits.lines[idx].second << " - " << seq2pos[w].second << " = "
//                              << cargo(findEdge(graph2, traits.lines[idx].second, seq2pos[w].second))
//                              << std::endl;
//                    ++nClosedLoops;
                    validInteractions[traits.lines[idx].first][seq2pos[w].first]
                        = cargo(findEdge(graph1, traits.lines[idx].first, seq2pos[w].first))
                        + cargo(findEdge(graph2, traits.lines[idx].second, seq2pos[w].second));

/*                    if (degree(graph1, traits.lines[j].first) > 0 && degree(graph2, traits.lines[j].second) > 0)
                    std::cout << "Interaction match: " << traits.lines[idx].first+1 << " - " << traits.lines[idx].second+1
                              << " | " << traits.lines[j].first+1 << " - " << traits.lines[j].second+1
                              << " // " << findEdge(graph1, traits.lines[idx].first, traits.lines[j].first)
                              << " - " << findEdge(graph2, traits.lines[idx].second, traits.lines[j].second)
                              << std::endl;
*/
                }
            }

        }
//        traits.numberOfSubgradients = nLoops + 1 - (nClosedLoops * 2); // TODO check this value: It should be referred to the alignment line only or to all the lambda?
    }
}

// ----------------------------------------------------------------------------
// Function computeBounds() version that make use of the lemon MWM
// ----------------------------------------------------------------------------
/*
void computeBounds(RnaAlignmentTraits & traits, InteractionScoreMap * validInteractions) // upper bound computation
{
    Graph<Undirected<double> > & graph1 = traits.bppGraphH.inter;
    Graph<Undirected<double> > & graph2 = traits.bppGraphV.inter;

    // Clear the weight of the upper bound
    for (std::size_t idx = 0; idx < length(traits.weightLineVect); ++idx)
    {
        traits.weightLineVect[idx].weight = 0; // reset best line score
    }

    // iterate all lines that are present in the alignment
    for (PositionPair const & line : traits.lines)
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
                if (validInteractions != NULL && line.first < value(adj_it1))
                {
                    for (PositionPair const & pairline : traits.lines)
                    { //TODO make this loop more efficient
                        if (pairline.first == value(adj_it1) && pairline.second == value(adj_it2))
                            (*validInteractions)[line.first][pairline.first] = edgeWeight;
                    }
                }

//                std::cerr << "Interaction match: " << line.first+1 << " - " << line.second+1
//                          << " | " << value(adj_it1)+1 << " - " << value(adj_it2)+1 << "\tprob = "
//                          << edgeWeight1 << " + " << edgeWeight-edgeWeight1 << "\n";

                // for upper bound do not care if interactions are closed by a line
                if (traits.weightLineVect[line.second].weight < edgeWeight / 2.0)
                {
                    traits.weightLineVect[line.second].weight = edgeWeight / 2.0;
                    traits.weightLineVect[line.second].seq1Index = line.first;
                    traits.weightLineVect[line.second].lineL.first = value(adj_it1);
                    traits.weightLineVect[line.second].lineL.second = value(adj_it2);
//                    std::cerr << "updated\n";
                }
            }
        }
    }
}
*/

double computeLowerBoundGreedy(InteractionScoreMap & interactions)
{
    // add vertices
    RnaInteractionGraph graph;
    forEach(interactions, [&graph] (ScMap const &) { addVertex(graph); });

    // add edges
    for (unsigned vertexIdx = 0; vertexIdx < length(interactions); ++vertexIdx)
        for (ScMap::iterator mapIt = interactions[vertexIdx].begin(); mapIt != interactions[vertexIdx].end(); ++mapIt)
            addEdge(graph, vertexIdx, mapIt->first, mapIt->second);

    return maximumWeightedMatchingGreedy<5>(graph);
};

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
// Function updateClosedLoops()
// ----------------------------------------------------------------------------

void updateClosedLoops(RnaAlignmentTraits & traits)
{
//    int tmpIndex = -1;
    for (PositionPair const & line : traits.lines)
    {
        if (traits.interactions[line.first].count(line.second) > 0)
        {
            RnaInteraction & interaction = traits.interactions[line.first][line.second];
            if (interaction.lineL.first ==
                traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.first
                && interaction.lineL.second ==
                   traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.second)
            {
                interaction.closedLoop = true;
/*                if (saveFoundInterPair) //TODO check if this lambda must be updated with the closed loop or not
                    // if the map is not empty means that the line have been already evaluated in a previous iteration
                {
                    traits.interactions[interaction.lineM.first][interaction.lineM.second].closedLoop = true;
                }
*/            }
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

    bool flagUpdateClosedLoops = false;
    // iterate all lines that are present in the alignment
    for (PositionPair const & line : traits.lines)
    {
//        std::cerr << "graph degree " << degree(graph1, line.first) << " - " << degree(graph2, line.second) << std::endl;
        // outgoing interactions in the first sequence
        bool proceed = false;
        if (degree(graph1, line.first) > 0 && degree(graph2, line.second) > 0)
        {
            if (traits.interactions[line.first].count(line.second) == 0)
            {
                proceed = true;
            } else if (traits.interactions[line.first][line.second].fromUBPairing)
            {
                proceed = true;
            }
            traits.interactions[line.first][line.second].iterUpdate = iter;
            if (proceed)
            {
                RnaInteraction & interaction = traits.interactions[line.first][line.second];
                findMaxWeight(interaction.maxProbScoreLine1, interaction.lineM.first, graph1, line.first);
                interaction.lineL.first = line.first;
//            std::cerr << interaction.lineL.first << " : " << interaction.lineM.first << " | " << interaction.maxProbScoreLine1
//                      << " - " << line.second << std::endl;
                findMaxWeight(interaction.maxProbScoreLine2, interaction.lineM.second, graph2, line.second);
                interaction.lineL.second = line.second;
//            std::cerr << interaction.lineL.second << " : " << interaction.lineM.second << " | " << interaction.maxProbScoreLine2
//                      << " - " << line.first << std::endl;
                interaction.weight = (interaction.maxProbScoreLine1 + interaction.maxProbScoreLine2) / 2;

//                std::cerr << interaction.lineL.first << " : " << interaction.lineL.second << " = "
//                          << interaction.weight << " | " << interaction.lineM.first << " : "
//                          << interaction.lineM.second << " (" << interaction.fromUBPairing << ")" << std::endl;
                interaction.fromUBPairing = false;
                flagUpdateClosedLoops = true;
                if (true) //TODO probably this must be the default and removed from the options
                    // if the map is not empty means that the line have been already evaluated in a previous iteration
                {
                    bool proceed2 = false;
                    if (traits.interactions[interaction.lineM.first].count(interaction.lineM.second) == 0)
                    {
                        proceed2 = true;
                    } else if (traits.interactions[interaction.lineM.first][interaction.lineM.second].weight < interaction.weight)
                    {
                        proceed2 = true;
                    }
                    //traits.interactions[interaction.lineM.first][interaction.lineM.second].iterUpdate = iter; //TODO check if needed
                    if (proceed2)
                    {
                        RnaInteraction & pair = traits.interactions[interaction.lineM.first][interaction.lineM.second];
                        pair.lineL.first = interaction.lineM.first;
                        pair.lineL.second = interaction.lineM.second;
                        pair.lineM.first = interaction.lineL.first;
                        pair.lineM.second = interaction.lineL.second;
                        pair.weight = interaction.weight;
                        pair.fromUBPairing = true;

//                        std::cerr << pair.lineL.first << " : " << pair.lineL.second << " = "
//                                  << pair.weight << " | " << pair.lineM.first << " : "
//                                  << pair.lineM.second << " (" << pair.fromUBPairing << ")" << std::endl;
                    }
                }
            }
        }
    }
    if (flagUpdateClosedLoops)
        updateClosedLoops(traits);
}

// ----------------------------------------------------------------------------
// Function computeNumberOfSubgradients()
// ----------------------------------------------------------------------------

void computeNumberOfSubgradients(RnaAlignmentTraits & traits, seqan::String<PositionPair> & unclosedLoops)
{
    traits.numberOfSubgradients = 0;
    for (PositionPair const & line : traits.lines)
    {
        if (traits.interactions[line.first].count(line.second) > 0)
        {
            bool stucturalLine = false;
/*
            RnaInteraction &lambdaL = traits.interactions[lineL.first][lineL.second];

            bool closedCircle = false;
            for (PositionPair const & lineM : traits.lines)
            {
                if (traits.interactions[lineM.first].count(lineM.second) > 0)
                    if (lambdaL.lineM.first == lineM.first && lambdaL.lineM.second == lineM.second)
                        closedCircle = true;
            }
            if (!closedCircle)
            {
                traits.numberOfSubgradients += 2;
                appendValue(unclosedLoops, lineL);
            }
*/
            RnaInteraction & interaction = traits.interactions[line.first][line.second];
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
                traits.numberOfSubgradients += 2;
                appendValue(unclosedLoops, line);
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function updateLambdaValues()
// ----------------------------------------------------------------------------

void updateLambdaValues(RnaAlignmentTraits & traits, seqan::String<PositionPair> const & unclosedLoops)
{
    for (PositionPair const & line : unclosedLoops)
    {
        if (traits.interactions[line.first].count(line.second) > 0)
        {
            RnaInteraction & interaction = traits.interactions[line.first][line.second];
            if (true || interaction.lineL.first !=
                traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.first
                || interaction.lineL.second !=
                   traits.interactions[interaction.lineM.first][interaction.lineM.second].lineM.second)
            {
                interaction.lambdaValue -= traits.stepSize;
                traits.interactions[interaction.lineM.first][interaction.lineM.second].lambdaValue += traits.stepSize;
            }
/*
            std::cout << line.first << " : " << line.second << " = " << interaction.weight << " -> "<< interaction.step << " | ";
            std::cout << interaction.lineM.first << " : "
                      << interaction.lineM.second << " = "
                      << traits.interactions[interaction.lineM.first][interaction.lineM.second].weight << " -> "
                      << traits.interactions[interaction.lineM.first][interaction.lineM.second].step
                      << std::endl;
*/
        }
    }
}

#endif //_INCLUDE_ALIGNMENT_EDGES_H_
