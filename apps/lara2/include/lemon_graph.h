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

void computeLowerBoundScore(TMapVect & lowerBound4Lemon, RnaAlignmentTraits & rnaAlign)
{
    lemon::SmartGraph lemonG;
    typedef lemon::SmartGraph::EdgeMap<TScoreValue > EdgeMap;
    typedef lemon::SmartGraph::NodeIt NodeIt;
    typedef lemon::SmartGraph::EdgeIt EdgeIt;
    EdgeMap weight(lemonG);
    // Add vertex for the LowerBoundGraph
    std::vector<lemon::SmartGraph::Node> fnv;
    for(unsigned i = 0; i < lowerBound4Lemon.size(); ++i)
    {
        lemon::SmartGraph::Node node = lemonG.addNode();
        fnv.push_back(node);
    }

    unsigned ll=0;
//    typedef Iterator<TLowerBoundGraph, AdjacencyIterator>::Type TAdjacencyIterator;
//    unsigned j;
    for(unsigned i = 0; i < lowerBound4Lemon.size(); ++i)
    { // With this function the direct graph storing all the couples of edges is generated
        for(auto const & trg_prob : lowerBound4Lemon[i])
        {
//        for (const auto& [trg, prob] : lowerBound4Lemon[i])
//            std::cout << "Planet " << name << ":\n" << description << "\n\n";
//        for (unsigned j = i+1; j < (lowerBound4Lemon[i].size()); ++j) {

//            std::cout << "(" << i << ":" << trg_prob.first << ") = " << trg_prob.second << std::endl;
            //lowerBound4Lemon[i][trg] << std::endl; //prob << std::endl;
            lemon::SmartGraph::Edge edge = lemonG.addEdge(fnv[i], fnv[trg_prob.first]);
//            lemon::SmartGraph::Edge edge = lemonG.addEdge(node1, node2);
//            weight[edge] = seqan::cargo(*it);
            weight[edge] = trg_prob.second; //lowerBound4Lemon[i][trg]; //prob;
            ++ll;
        }
    }

    std::cerr << "Number of edges = " << ll << std::endl;
    std::cerr << "Nodes:";
    for (NodeIt i(lemonG); i!=lemon::INVALID; ++i)
        std::cerr << " " << lemonG.id(i);
    std::cerr << std::endl;

    std::cerr << "Edges:";
    for (EdgeIt i(lemonG); i!=lemon::INVALID; ++i)
        std::cerr << " (" << lemonG.id(lemonG.u(i)) << "," << lemonG.id(lemonG.v(i)) << ")";
    std::cerr << std::endl;

    // Do stuff
    lemon::MaxWeightedMatching<lemon::SmartGraph, EdgeMap> mwm(lemonG, weight);
    mwm.run();

    rnaAlign.lowerLemonBound.mwmPrimal = mwm.matchingWeight();
    rnaAlign.lowerLemonBound.mwmDual = mwm.dualValue();
    rnaAlign.lowerLemonBound.mwmCardinality = mwm.matchingSize();
    std::cerr << "Lemon result: weight= " << mwm.matchingWeight() << "\tdual= " << mwm.dualValue()
              << "\tsize= " << mwm.matchingSize() << "\n";


    for(lemon::SmartGraph::EdgeIt e(lemonG); e!=lemon::INVALID; ++e){
        if(mwm.matching(e)){
            std::cerr << "weight(" << lemonG.id(lemonG.u(e)) << ","
                      << lemonG.id(lemonG.v(e)) << ")="<<weight[e]<<std::endl;
        }
    }

//    createLemonGraph(options, rnaAlign, lemonG);
};

}

#endif //_INCLUDE_LEMON_GRAPH_H_
