// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Gianvito Urgese <gianvito.urgese@gmail.com>
// ==========================================================================

#ifndef SEQAN_INCLUDE_ALIGN_RNA_SCORE_STRUCTURE_RNA_H_
#define SEQAN_INCLUDE_ALIGN_RNA_SCORE_STRUCTURE_RNA_H_

#include "data_types.h"

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

namespace seqan
{

template<typename TScoreMatrix, typename TOutgoingInteractions>
struct RnaStructureScore;

// ----------------------------------------------------------------------------
// Class RnaStructureScore Score
// ----------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions>  //typename TGap
class Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> >
{
public:
    TScoreMatrix matrix;
    TOutgoingInteractions *interactions;

    TValue getMapLineValue(unsigned seq1_pos, unsigned seq2_pos) const
    {
        if ((*interactions)[seq1_pos].count(seq2_pos) > 0)
        {
//            if((*interactions)[seq1_pos][seq2_pos].closedLoop)
//                return (*interactions)[seq1_pos][seq2_pos].weight;
//            else
            return ((*interactions)[seq1_pos][seq2_pos].lambdaValue + (*interactions)[seq1_pos][seq2_pos].weight);
        }
        else
        {
            return 0;
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                     [PositionSeq Score]
// --------------------------------------------------------------------------

// Returns the type that holds a sequence entry.  This is used for abstracting away the access to sequence characters.
template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSequence>
struct SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> >, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSequence>
struct SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const, TSequence> :
    SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> >, TSequence>
{
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                      [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TScoreValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const & /*sScheme*/,
                      TSequence const &seq, TPosition pos)
{
//  std::cout <<"ConsensusScoreSequenceEntry<TSequence>(seq, pos) "<< ConsensusScoreSequenceEntry<TSequence>(seq, pos)
//  std::cout << " seq " << seq << " pos " << pos << " sequenceEntryForScore " << std::endl;
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

// --------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()                   [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const &entry1,
    ConsensusScoreSequenceEntry<TSeq2> const &entry2)
{
//  unsigned posit=position(entry1);
//  if ((int)position(entry1) == 0)
//  std::cout << "entry1._seq[entry1._pos] = "<< (*entry1._seq)[posit] << " entry1._pos = " << posit << "  ||  ";
//  std::cout << "entry2._seq[entry2._pos] = "<< (*entry2._seq)[position(entry2)] << " entry2._pos = " << position(entry2) << " scoreGapExtendHorizontal " << std::endl;
//  std::cout << (int)value(entry1) << " " << (*entry1._seq)[posit] << std::endl;
    return scoreGapExtendHorizontal(me.matrix, (*entry1._seq)[position(entry1)],
                                    (*entry2._seq)[position(entry2)]);
//    return scoreGapExtendHorizontal(me.score_matrix, (unsigned)ordValue(entry1._seq[0][entry1._pos]), (unsigned)ordValue(entry2._seq[0][entry2._pos]));
//    return scoreGapExtendHorizontal(me.score_matrix, Nothing(), Nothing());
}

// --------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()                     [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const &entry1,
    ConsensusScoreSequenceEntry<TSeq2> const &entry2)
{
//  return scoreGapOpenHorizontal(me.score_matrix, (unsigned)ordValue(entry1._seq[0][entry1._pos]), (unsigned)ordValue(entry2._seq[0][entry2._pos]));
    return scoreGapOpenHorizontal(me.matrix, (*entry1._seq)[position(entry1)],
                                  (*entry2._seq)[position(entry2)]);
}

// --------------------------------------------------------------------------
// Function scoreGapExtendVertical()                     [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
// return scoreGapOpenVertical(me.score_matrix, entry1._seq[0][entry1.pos], entry2._seq[0][entry2.pos]);
    return scoreGapOpenVertical(me.matrix, Nothing(), Nothing());
}

// --------------------------------------------------------------------------
// Function scoreGapOpenVertical()                       [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const &me,
    ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
    ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
// return scoreGapExtendVertical(me.score_matrix, entry1._seq[0][entry1.pos], entry2._seq[0][entry2.pos]);
    return scoreGapExtendVertical(me.matrix, Nothing(), Nothing());
}



// --------------------------------------------------------------------------
// Function score()                                      [RnaStructure Score]
// --------------------------------------------------------------------------

template<typename TValue, typename TScoreMatrix, typename TOutgoingInteractions, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, RnaStructureScore<TScoreMatrix, TOutgoingInteractions> > const &me,
      ConsensusScoreSequenceEntry<TSeq1> const &entry1,
      ConsensusScoreSequenceEntry<TSeq2> const &entry2)
{
    /*
    if (me.getMapLineValue(position(entry1),position(entry2)) != 0)
    {
        std::cout << (*entry1._seq)[position(entry1)] << " " << (*entry2._seq)[position(entry2)] << " "
                  << position(entry1) << " " << position(entry2) << " "
                  << me.getMapLineValue(position(entry1), position(entry2)) << " "
                  << score(me.score_matrix, (*entry1._seq)[position(entry1)], (*entry2._seq)[position(entry2)])
                  << std::endl;
    }
    */
    // " mapLine =  " << me._mapLine[position(entry1)][position(entry2)] << std::endl; // me._mapLine[position(entry1)][position(entry2)]
// return score(me.score_matrix, (*entry1._seq)[position(entry1)] , (*entry2._seq)[position(entry2)]); // Normal Score using the substitutional matrix
    return score(me.matrix, (*entry1._seq)[position(entry1)],
                 (*entry2._seq)[position(entry2)]) + me.getMapLineValue(position(entry1), position(entry2));
}

} // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_ALIGN_RNA_SCORE_STRUCTURE_RNA_H_
