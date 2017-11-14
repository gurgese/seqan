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
// This file contains the
// ==========================================================================

#ifndef _INCLUDE_LEMON_GRAPH_H_
#define _INCLUDE_LEMON_GRAPH_H_

// ----------------------------------------------------------------------------
// Lemon headers
// ----------------------------------------------------------------------------

#include <lemon/list_graph.h>
#include <lemon/smart_graph.h>
#include <lemon/matching.h>

// ----------------------------------------------------------------------------
// Lara data types
// ----------------------------------------------------------------------------

#include "data_types.h"

// ============================================================================
// Functions
// ============================================================================

namespace myLemon {

// ----------------------------------------------------------------------------
// Function computeLowerBound()
// ----------------------------------------------------------------------------

double computeLowerBoundScore(InteractionScoreMap const & validInteractionScores)
{
    typedef lemon::SmartGraph::EdgeMap<double> EdgeMap;
    lemon::SmartGraph lemonG;

    // Add nodes for the LowerBoundGraph
    std::vector<lemon::SmartGraph::Node> nodes;
    forEach(validInteractionScores, [&lemonG, &nodes] (ScMap const &) { nodes.push_back(lemonG.addNode()); });

    // Add the edges and assign weights.
    EdgeMap weight(lemonG);
    for (unsigned nodeIdx = 0; nodeIdx < validInteractionScores.size(); ++nodeIdx)
    {
        for (std::pair<unsigned, double> const & pairedProb : validInteractionScores[nodeIdx])
            weight[lemonG.addEdge(nodes[nodeIdx], nodes[pairedProb.first])] = pairedProb.second;
    }

    // Perform the MWM.
    lemon::MaxWeightedMatching<lemon::SmartGraph, EdgeMap> mwm(lemonG, weight);
    mwm.run();
    return mwm.matchingWeight();
};

}

#endif //_INCLUDE_LEMON_GRAPH_H_
