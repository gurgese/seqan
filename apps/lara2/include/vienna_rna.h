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

#ifndef _INCLUDE_VIENNA_RNA_H_
#define _INCLUDE_VIENNA_RNA_H_

// ----------------------------------------------------------------------------
// Vienna headers
// ----------------------------------------------------------------------------

extern "C" {
    #include <ViennaRNA/data_structures.h>
    #include <ViennaRNA/params.h>
    #include <ViennaRNA/utils.h>
    #include <ViennaRNA/eval.h>
    #include <ViennaRNA/fold.h>
    #include <ViennaRNA/part_func.h>
    #include <ViennaRNA/PS_dot.h>
}

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <omp.h>

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeBppMatrix()
// ----------------------------------------------------------------------------

void computeBppMatrix(RnaRecord & rnaRecord, double thrBppm, bool logStructureScoring)
{
    RnaStructureGraph bppMatrGraph;

    // apply default model details
    vrna_md_t md_p{};
    set_model_details(&md_p);

    // get a vrna_fold_compound with MFE and PF DP matrices and default model details
    String<char, CStyle> seq = rnaRecord.sequence;
    vrna_fold_compound_t *vc = vrna_fold_compound(toCString(seq), &md_p, VRNA_OPTION_MFE | VRNA_OPTION_PF);

    // save gibbs free energy
    std::vector<char> structure;
    structure.resize(length(rnaRecord.sequence) + 1);
    float gibbs = vrna_pf(vc, structure.data());
    bppMatrGraph.energy = gibbs;

    vrna_plist_t *pl1;
    pl1 = vrna_plist_from_probs(vc, thrBppm);

    // get size of pl1
    unsigned size;
    vrna_plist_t *ptr = pl1;
    for (size = 0u; ptr->i; ++size, ++ptr);

    // define string for graph specs
    //TODO this data must be formatted in a smart way
    bppMatrGraph.specs = "vrna_fold_compound(<Sequence>, <Vienna Model Details>, VRNA_OPTION_MFE | VRNA_OPTION_PF)";

    // add vertices to graph
    for (unsigned idx=0; idx < length(rnaRecord.sequence); ++idx)
        addVertex(bppMatrGraph.inter);

    if (logStructureScoring)
    {
        /*
        double minProb = 1.0;
        for (unsigned i = 0; i < size; ++i)
        {
            if (pl1[i].p < minProb)
                minProb = pl1[i].p;
        }
         */
        double const minProb = 0.003; // taken from LISA > Lara
        for (unsigned i = 0; i < size; ++i)
        {
            SEQAN_ASSERT(pl1[i].i > 0 && static_cast<unsigned>(pl1[i].i) <= length(rnaRecord.sequence));
            SEQAN_ASSERT(pl1[i].j > 0 && static_cast<unsigned>(pl1[i].j) <= length(rnaRecord.sequence));
            // convert indices from range 1..length to 0..length-1
            if (pl1[i].p > minProb)
                addEdge(bppMatrGraph.inter, pl1[i].i - 1, pl1[i].j - 1, log(pl1[i].p/minProb));
        }
    }
    else
    {
        for (unsigned i = 0; i < size; ++i)
        {
            SEQAN_ASSERT(pl1[i].i > 0 && static_cast<unsigned>(pl1[i].i) <= length(rnaRecord.sequence));
            SEQAN_ASSERT(pl1[i].j > 0 && static_cast<unsigned>(pl1[i].j) <= length(rnaRecord.sequence));
            // convert indices from range 1..length to 0..length-1
            addEdge(bppMatrGraph.inter, pl1[i].i - 1, pl1[i].j - 1, pl1[i].p);
        }
    }
    append(rnaRecord.bppMatrGraphs, bppMatrGraph);

    /*
    // fixed graph: add specs and energy
    RnaStructureGraph fixedGraph;
    fixedGraph.specs = "vrna_fold_compound(<Sequence>, <Vienna Model Details>, VRNA_OPTION_MFE | VRNA_OPTION_PF)";
    fixedGraph.energy = gibbs;
    append(rnaRecord.fixedGraphs, fixedGraph);
    */

    // free memory occupied by vrna_fold_compound
    vrna_fold_compound_free(vc);
}

#endif //_INCLUDE_VIENNA_RNA_H_
