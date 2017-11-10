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

void computeMissingInteractions(TRnaVect & rnaRecordVector, LaraOptions const & options)
{
    bool const logStructureScoring = options.structureScoring == LOGARITHMIC;

//#pragma omp parallel for num_threads(options.threads)
    for (typename Size<TRnaVect>::Type idx = 0; idx < length(rnaRecordVector); ++idx)
    {
        if (empty(rnaRecordVector[idx].bppMatrGraphs))  // if dotplot or extended bpseq data are not present
            computeBppMatrix(rnaRecordVector[idx], options.thrBppm, logStructureScoring);
    }
}

// ----------------------------------------------------------------------------
// Function setScoreMatrix()
// ----------------------------------------------------------------------------

void setScoreMatrix(LaraOptions & options)
{
    options.laraScoreMatrix.data_gap_extend = options.laraGapExtend;
    options.laraScoreMatrix.data_gap_open   = options.laraGapOpen;
    if (empty(options.laraScoreMatrixName))
    {
        _VV(options, "Predefined RIBOSUM matrix will be used");
        setDefaultScoreMatrix(options.laraScoreMatrix, TRibosum());
    }
    else if (loadScoreMatrix(options.laraScoreMatrix, toCString(options.laraScoreMatrixName)))
    {
        _VV(options, "Provided scoring matrix will be used " << options.laraScoreMatrixName);
    }
    else
    {
        std::cerr << "Matrix file could not be opened: " << options.laraScoreMatrixName
                  << "). Predefined RIBOSUM matrix will be used." << std::endl;
        setDefaultScoreMatrix(options.laraScoreMatrix, TRibosum());
    }

    // scale the matrix
    options.laraScoreMatrix.data_gap_extend /= options.sequenceScale;
    options.laraScoreMatrix.data_gap_open   /= options.sequenceScale;
    for (double & matrixEntry : options.laraScoreMatrix.data_tab)
        matrixEntry /= options.sequenceScale;
}

// ----------------------------------------------------------------------------
// Function firstSimdAlignsGlobalLocal()
// ----------------------------------------------------------------------------

void initialAlignment(String<double> & results, StringSet<RnaAlignment> & alignments, LaraOptions const & options)
{
    TScoreMatrix scoreMatrix = options.laraScoreMatrix;
    scoreMatrix.data_gap_extend = options.generatorGapExtend / options.sequenceScale;
    scoreMatrix.data_gap_open   = options.generatorGapOpen   / options.sequenceScale;

    if (!options.globalLocal)  //TODO implement the global-unconstrained alignment using the parameters in the options
    {
        if (options.affineLinearDgs == 0)
            results = globalAlignment(alignments, scoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            results = globalAlignment(alignments, scoreMatrix, LinearGaps());
        else
            results = globalAlignment(alignments, scoreMatrix, DynamicGaps());

    }
    else
    {
        if (options.affineLinearDgs == 0)
            results = localAlignment(alignments, scoreMatrix, AffineGaps());
        else if (options.affineLinearDgs == 1)
            results = localAlignment(alignments, scoreMatrix, LinearGaps());
        else
            results = localAlignment(alignments, scoreMatrix, DynamicGaps());
    }
};

// ----------------------------------------------------------------------------
// Function simdAlignsGlobalLocal()
// ----------------------------------------------------------------------------
//TODO modify this function in order to implement the SIMD alignment
void structuralAlignment(String<double> & results, StringSet<RnaAlignment> & alignments,
                         RnaAlignmentTraitsVector const & alignmentTraits, LaraOptions const & options)
{
    if (!options.globalLocal)  //TODO implement the global-unconstrained alignment using the parameters in the options
    {
        if (options.affineLinearDgs == 0)
        {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = globalAlignment(alignments[idx], alignmentTraits[idx].structureScore, AffineGaps());

        }
        else if (options.affineLinearDgs == 1)
        {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = globalAlignment(alignments[idx], alignmentTraits[idx].structureScore, LinearGaps());

        }
        else {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = globalAlignment(alignments[idx], alignmentTraits[idx].structureScore, DynamicGaps());

        }
    }
    else
    {
        if (options.affineLinearDgs == 0)
        {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = localAlignment(alignments[idx], alignmentTraits[idx].structureScore, AffineGaps());

        }
        else if (options.affineLinearDgs == 1)
        {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = localAlignment(alignments[idx], alignmentTraits[idx].structureScore, LinearGaps());

        }
        else {
//#pragma omp parallel for num_threads(options.threads)
            for (unsigned idx = 0; idx < length(alignments); ++idx) // TODO replace this function with the SIMD implementation for execute in PARALLEL
                results[idx] = localAlignment(alignments[idx], alignmentTraits[idx].structureScore, DynamicGaps());

        }
    }
};

// ----------------------------------------------------------------------------
// Function crossproduct()
// ----------------------------------------------------------------------------

inline void _fillVectors (RnaSeqSet & setH, RnaSeqSet & setV, RnaAlignmentTraitsVector::iterator & alignInfo,
                          TRnaVect::iterator const & it1, TRnaVect::iterator const & it2)
{
    appendValue(setH, it1->sequence);
    appendValue(setV, it2->sequence);
    alignInfo->bppGraphH = front(it1->bppMatrGraphs);
    alignInfo->bppGraphV = front(it2->bppMatrGraphs);
    SEQAN_ASSERT_EQ(length(back(setH)), numVertices(alignInfo->bppGraphH.inter));
    SEQAN_ASSERT_EQ(length(back(setV)), numVertices(alignInfo->bppGraphV.inter));
    ++alignInfo;
}

// unique combination of all sequences of one single set
void crossproduct(RnaSeqSet & setH, RnaSeqSet & setV, RnaAlignmentTraitsVector & alignmentTraits, TRnaVect & seqs)
{
    typename Size<TRnaVect>::Type const len = length(seqs);
    if (len == 0)
        return;

    reserve(setH, len * (len - 1) / 2);
    reserve(setV, len * (len - 1) / 2);
    resize(alignmentTraits, len * (len - 1) / 2);
    RnaAlignmentTraitsVector::iterator alignInfo = alignmentTraits.begin();

    unsigned p = 0;
    for (TRnaVect::iterator it1 = seqs.begin(); it1 != seqs.end(); ++it1)
    {
        for (TRnaVect::iterator it2 = it1 + 1u; it2 != seqs.end(); ++it2)
        {
            _fillVectors(setH, setV, alignInfo, it1, it2);
            alignmentTraits[p].idBppSeqH = static_cast<unsigned>(std::distance(seqs.begin(), it1));
            alignmentTraits[p].idBppSeqV = static_cast<unsigned>(std::distance(seqs.begin(), it2));
            ++p;
        }
    }
}

// combination of all sequences of two sets
void crossproduct(RnaSeqSet & setH, RnaSeqSet & setV, RnaAlignmentTraitsVector & alignmentTraits,
                  TRnaVect & seqs1, TRnaVect & seqs2) {
    if (empty(seqs2))
    {
        crossproduct(setH, setV, alignmentTraits, seqs1);
    }
    else
    {
        reserve(setH, length(seqs1) * length(seqs2));
        reserve(setV, length(seqs1) * length(seqs2));
        resize(alignmentTraits, length(seqs1) * length(seqs2));
        RnaAlignmentTraitsVector::iterator alignInfo = alignmentTraits.begin();

        unsigned p = 0;
        for (TRnaVect::iterator it1 = seqs1.begin(); it1 != seqs1.end(); ++it1)
        {
            for (TRnaVect::iterator it2 = seqs2.begin(); it2 != seqs2.end(); ++it2)
            {
                _fillVectors(setH, setV, alignInfo, it1, it2);
                alignmentTraits[p].idBppSeqH = std::distance(seqs1.begin(), it1);
                alignmentTraits[p].idBppSeqV = std::distance(seqs2.begin(), it2);
                ++p;
            }
        }
    }
}

void prepareAlignmentVector(StringSet<RnaAlignment> & alignments, RnaAlignmentTraitsVector & alignmentTraits,
                            RnaStructContentsPair & filecontents)
{
    RnaSeqSet setH;
    RnaSeqSet setV;
    crossproduct(setH, setV, alignmentTraits, filecontents.first.records, filecontents.second.records);
    SEQAN_ASSERT_EQ(length(setH), length(setV));
    SEQAN_ASSERT_EQ(length(setH), length(alignmentTraits));

    resize(alignments, length(setH));
    auto it = makeZipIterator(begin(setH), begin(setV), begin(alignments));
    auto itEnd = makeZipIterator(end(setH), end(setV), end(alignments));

    while (it != itEnd)
    {
        RnaAlignment align;
        resize(rows(align), 2);
        assignSource(row(align, 0), std::get<0>(*it));
        assignSource(row(align, 1), std::get<1>(*it));
        std::get<2>(*it) = align;
        ++it;
    }
    SEQAN_ASSERT_EQ(length(alignments), length(alignmentTraits));
}

#endif //_INCLUDE_STRUCT_ALIGN_H_
