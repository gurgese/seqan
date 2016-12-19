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
// Authors: Gianvito Urgese <gianvito.urgese@polito.it>
//          Joerg Winkler <j.winkler@fu-berlin.de>
// ==========================================================================
// This file contains
// ==========================================================================
#ifndef _INCLUDE_STRUCT_ALIGN_H_
#define _INCLUDE_STRUCT_ALIGN_H_

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

#include "vienna_rna.h"

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function bppInteractionGraphBuild()
// ----------------------------------------------------------------------------

template <typename TOption>
void bppInteractionGraphBuild(TRnaVect & rnaSeqs, TOption const & options)
{
#pragma omp parallel for num_threads(options.threads)
    for (typename Size<TRnaVect>::Type i = 0; i < length(rnaSeqs); ++i)
    {
        if (empty(rnaSeqs[i].bppMatrGraphs))  // if dotplot or extended bpseq data are not present
            computeBppMatrix(rnaSeqs[i], options);
    }
}

// ----------------------------------------------------------------------------
// Function setScoreMatrix()
// ----------------------------------------------------------------------------

template <typename TOptions>
void setScoreMatrix(TOptions & options)
{
    options.laraScoreMatrix.data_gap_extend = options.laraGapExtend;
    options.laraScoreMatrix.data_gap_open = options.laraGapOpen;
    if (options.laraScoreMatrixName != "")
    {
        loadScoreMatrix(options.laraScoreMatrix, toCString(getAbsolutePath(toCString(options.laraScoreMatrixName))));
        _V(options, "Provided scoring matrix will be used " << options.laraScoreMatrixName);
//        showScoringMatrix(options.laraScoreMatrix);
    }
    else
    {
        _V(options, "Predefined RIBOSUM matrix will be used");
        setDefaultScoreMatrix(options.laraScoreMatrix, TRibosum());
//        showScoringMatrix(options.laraScoreMatrix);
    }
}

// ----------------------------------------------------------------------------
// Function firstSimdAlignsGlobalLocal()
// ----------------------------------------------------------------------------

template <typename TResultsSimd, typename TAlignsSimd, typename TOptions>
void firstSimdAlignsGlobalLocal(TResultsSimd & resultsSimd, TAlignsSimd & alignsSimd, TOptions const & options)
{
    if (!options.globalLocal)  //TODO implement the global-unconstrained alignment using the parameters in the options
    {
        if (options.affineLinearDgs == 0)
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, LinearGaps());
        else
            resultsSimd = globalAlignment(alignsSimd, options.laraScoreMatrix, DynamicGaps());

    } else
    {
        if (options.affineLinearDgs == 0)
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, LinearGaps());
        else
            resultsSimd = localAlignment(alignsSimd, options.laraScoreMatrix, DynamicGaps());
    }
};

template <typename TAlignsSimd>
void createSimdAligns(TAlignsSimd & alignsSimd, RnaSeqSet const & setH, RnaSeqSet const & setV)
{
    resize(alignsSimd, length(setH));
    auto it = makeZipIterator(begin(setH), begin(setV), begin(alignsSimd));
    auto itEnd = makeZipIterator(end(setH), end(setV), end(alignsSimd));

    while (it != itEnd)
    {
        TAlign align;
        resize(rows(align), 2);
        assignSource(row(align, 0), std::get<0>(*it));
        assignSource(row(align, 1), std::get<1>(*it));
        std::get<2>(*it) = align;
        ++it;
    }
}

// ----------------------------------------------------------------------------
// Function crossproduct()
// ----------------------------------------------------------------------------

inline void _fillVectors (RnaSeqSet & setH, RnaSeqSet & setV, TRnaAlignVect::iterator & alignInfo,
                          TRnaVect::iterator const & it1, TRnaVect::iterator const & it2)
{
    // in this way the alignment map structure will be always created with the maximum size
    if (length(it1->sequence) > length(it2->sequence))
    {
        appendValue(setH, it1->sequence);
        appendValue(setV, it2->sequence);
        alignInfo->bppGraphH = it1->bppMatrGraphs[alignInfo->idBppSeqH];
        alignInfo->bppGraphV = it2->bppMatrGraphs[alignInfo->idBppSeqV];
    }
    else
    {
        appendValue(setH, it2->sequence);
        appendValue(setV, it1->sequence);
        alignInfo->bppGraphH = it2->bppMatrGraphs[alignInfo->idBppSeqH];
        alignInfo->bppGraphV = it1->bppMatrGraphs[alignInfo->idBppSeqV];
    }
    ++alignInfo;
}

// unique combination of all sequences of one single set
void crossproduct(RnaSeqSet & setH, RnaSeqSet & setV, TRnaAlignVect & rnaAligns, TRnaVect & seqs)
{
    typename Size<TRnaVect>::Type const len = length(seqs);
    if (len == 0)
        return;

    reserve(setH, len * (len - 1) / 2);
    reserve(setV, len * (len - 1) / 2);
    resize(rnaAligns, len * (len - 1) / 2);
    TRnaAlignVect::iterator alignInfo = rnaAligns.begin();

    for (TRnaVect::iterator it1 = seqs.begin(); it1 != seqs.end(); ++it1)
        for (TRnaVect::iterator it2 = it1 + 1u; it2 != seqs.end(); ++it2)
            _fillVectors(setH, setV, alignInfo, it1, it2);
}

// combination of all sequences of two sets
void crossproduct(RnaSeqSet & setH, RnaSeqSet & setV, TRnaAlignVect & rnaAligns,
                  TRnaVect & seqs1, TRnaVect & seqs2)
{
    if (empty(seqs2))
    {
        crossproduct(setH, setV, rnaAligns, seqs1);
        return;
    }

    reserve(setH, length(seqs1) * length(seqs2));
    reserve(setV, length(seqs1) * length(seqs2));
    resize(rnaAligns, length(seqs1) * length(seqs2));
    TRnaAlignVect::iterator alignInfo = rnaAligns.begin();

    for (TRnaVect::iterator it1 = seqs1.begin(); it1 != seqs1.end(); ++it1)
        for (TRnaVect::iterator it2 = seqs2.begin(); it2 != seqs2.end(); ++it2)
            _fillVectors(setH, setV, alignInfo, it1, it2);
}

#endif //_INCLUDE_STRUCT_ALIGN_H_