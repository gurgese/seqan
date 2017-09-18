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

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TScoreMatrix, typename TLambVect>
struct RnaStructureScore;

// ----------------------------------------------------------------------------
// Class PositionSeqScore Score
// ----------------------------------------------------------------------------

//template <typename TValue, typename TMapline, typename TSequence, typename TAlign>
//class Score<TValue, PositionSeqScore>
template <typename TValue, typename TScoreMatrix, typename TLambVect>  //typename TGap
class Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> >
//template <typename TValue, typename TSequence>  //typename TGap
//class Score<TValue, RnaStructureScore<TSequence> >
{
public:

//	typedef typename Value<TSequence>::Type TSeqValue;
    // A table of position x (ord value) giving the counts of the characters at the given positions.
//    String<std::map<TPosition, TValue> > mapline;
//	TScoreMatrix _scoreMatrix;
// Static computation of the required array size.
//enum {
//	VALUE_SIZE = ValueSize<typename Value<TSequence>::Type>::VALUE,
//	TAB_SIZE = VALUE_SIZE * VALUE_SIZE
//};
    // The data table.
//    Score<TValue, ScoreMatrix<TSeqValue, Ribosume65> >
    TScoreMatrix score_matrix;
//    TValue _score_matrix_tab[TAB_SIZE];

//	Score<TValue, ScoreMatrix<Rna, Default> > scoringSchemeMatrix; //(data_gap_extend, data_gap_open)
//	String<std::map<size_t, TValue> >   _mapLine; //Rene
//	TMapline   _mapLine;
    String<std::map<unsigned, TValue> > * _mapLine;
    TLambVect * lamb;

    TValue getMapLineValue(unsigned seq1_pos, unsigned seq2_pos) const
    {
 //   	if (((*_mapLine)[seq1_pos]).find(seq2_pos) !=  ((*_mapLine)[seq1_pos]).end())
		if ((*lamb)[seq1_pos].map.count(seq2_pos) == 1)
    		return ((*lamb)[seq1_pos].map[seq2_pos].step + (*lamb)[seq1_pos].map[seq2_pos].maxProbScoreLine);
    	else
    		return 0;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction SequenceEntryForScore                     [PositionSeq Score]
// --------------------------------------------------------------------------

// Returns the type that holds a sequence entry.  This is used for abstracting away the access to sequence characters.
template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSequence>
struct SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> >, TSequence>
{
    typedef ConsensusScoreSequenceEntry<TSequence> Type;
};

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSequence>
struct SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> > const, TSequence> :
            SequenceEntryForScore<Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> >, TSequence>
{};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function fillUpdateMapline()                          [RnaStructure Score]
// --------------------------------------------------------------------------
//TODO to be tested deeply
/*template <typename TAlign, typename TMapLine, typename TSequence, typename TValue>
void fillUpdateMapline(TMapLine & mapline, TAlign const & align, TSequence const & seq1, TValue const & factor)
//FIXME the sequence is used to pass the sequence type only
{
	std::cerr << "function fillUpdateMapline\n";
	TValue struct_score = 0.1 * factor; // this score must be updated including structural score

    seqan::Gaps<TSequence> row0 = row(align, 0);
    seqan::Gaps<TSequence> row1 = row(align, 1);
//    std::cout << "This is " << row1 << std::endl;
    unsigned row0Begin = clippedBeginPosition(row(align, 0));
    unsigned row1Begin = clippedBeginPosition(row(align, 1));
    unsigned gap0 = 0;
    unsigned gap1 = 0;
    std::cout << "row0Begin " << row0Begin << " row1Begin " << row1Begin << std::endl;
    for (unsigned i = 0; i < length(row(align, 0)); ++i)  //maximum size of this string is length(mapline)
	{
		if (row0[i] == '-')
		{
			++gap0;
		}
        else if (row1[i] == '-')
        {
			++gap1;
		}
		if(row0[i] == row1[i])
		{
			mapline[i + row0Begin - gap0][i + row1Begin - gap1] = struct_score * i;
            // must be decided if the computed score must be summed or just replaced
			std::cout << i + row0Begin - gap0 << ":" << i + row1Begin - gap1 << "/"
                      << mapline[i + row0Begin - gap0][i + row1Begin - gap1] << "\t";
		}
	}
    std::cout << "Sequence 1 " << row(align, 1) << std::endl;
//    myMapType::const_iterator it=myMap.begin(); it!=myMap.end(); ++it
    typedef typename Value<TMapLine>::Type TMap;
    for (unsigned i = 0; i < seqan::length(mapline); ++i)  //maximum size of this string is length(mapline)
	{
    	for (typename TMap::const_iterator it = mapline[i].begin(); it != mapline[i].end(); ++it)
        { // in this way only full lists are printed
            std::cout << i << ":" << it->first << " " << "/" << it->second << "\n";
        }
	}
}
*/
//template <typename TScore>
//void showScoringMatrix(TScore sSchemeRna)
//{
//	for(unsigned i=0; i< length(sSchemeRna.score_matrix_tab); ++i)
//		std::cout << sSchemeRna.score_matrix_tab[i] << "\t";
//	std::cout << "\n" << length(sSchemeRna.score_matrix_tab) << std::endl;
//
//}

// --------------------------------------------------------------------------
// Function sequenceEntryForScore()                      [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TScoreValue, typename TScoreMatrix, typename TLambVect, typename TSequence, typename TPosition>
inline ConsensusScoreSequenceEntry<TSequence>
sequenceEntryForScore(Score<TScoreValue, RnaStructureScore<TScoreMatrix, TLambVect> > const & /*sScheme*/,
                      TSequence const & seq, TPosition pos)
{
//	std::cout <<"ConsensusScoreSequenceEntry<TSequence>(seq, pos) "<< ConsensusScoreSequenceEntry<TSequence>(seq, pos)
//	std::cout << " seq " << seq << " pos " << pos << " sequenceEntryForScore " << std::endl;
    return ConsensusScoreSequenceEntry<TSequence>(seq, pos);
}

// --------------------------------------------------------------------------
// Function scoreGapExtendHorizontal()                   [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
        Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> > const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
//	unsigned posit=position(entry1);
//	if ((int)position(entry1) == 0)
//	std::cout << "entry1._seq[entry1._pos] = "<< (*entry1._seq)[posit] << " entry1._pos = " << posit << "  ||  ";
//	std::cout << "entry2._seq[entry2._pos] = "<< (*entry2._seq)[position(entry2)] << " entry2._pos = " << position(entry2) << " scoreGapExtendHorizontal " << std::endl;
//	std::cout << (int)value(entry1) << " " << (*entry1._seq)[posit] << std::endl;
	return scoreGapExtendHorizontal(me.score_matrix, (*entry1._seq)[position(entry1)] ,
                                    (*entry2._seq)[position(entry2)]);
//    return scoreGapExtendHorizontal(me.score_matrix, (unsigned)ordValue(entry1._seq[0][entry1._pos]), (unsigned)ordValue(entry2._seq[0][entry2._pos]));
//    return scoreGapExtendHorizontal(me.score_matrix, Nothing(), Nothing());
}

// --------------------------------------------------------------------------
// Function scoreGapOpenHorizontal()                     [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
        Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> > const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & entry1,
        ConsensusScoreSequenceEntry<TSeq2> const & entry2)
{
//	return scoreGapOpenHorizontal(me.score_matrix, (unsigned)ordValue(entry1._seq[0][entry1._pos]), (unsigned)ordValue(entry2._seq[0][entry2._pos]));
	return scoreGapOpenHorizontal(me.score_matrix, (*entry1._seq)[position(entry1)] ,
                                  (*entry2._seq)[position(entry2)]);
}

// --------------------------------------------------------------------------
// Function scoreGapExtendVertical()                     [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
        Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> > const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
//	return scoreGapOpenVertical(me.score_matrix, entry1._seq[0][entry1.pos], entry2._seq[0][entry2.pos]);
	return scoreGapOpenVertical(me.score_matrix, Nothing(), Nothing());
}

// --------------------------------------------------------------------------
// Function scoreGapOpenVertical()                       [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
        Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> > const & me,
        ConsensusScoreSequenceEntry<TSeq1> const & /*entry1*/,
        ConsensusScoreSequenceEntry<TSeq2> const & /*entry2*/)
{
//	return scoreGapExtendVertical(me.score_matrix, entry1._seq[0][entry1.pos], entry2._seq[0][entry2.pos]);
	return scoreGapExtendVertical(me.score_matrix, Nothing(), Nothing());
}



// --------------------------------------------------------------------------
// Function score()                                      [RnaStructure Score]
// --------------------------------------------------------------------------

template <typename TValue, typename TScoreMatrix, typename TLambVect, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, RnaStructureScore<TScoreMatrix, TLambVect> >  const & me,
      ConsensusScoreSequenceEntry<TSeq1> const & entry1,
      ConsensusScoreSequenceEntry<TSeq2> const & entry2)
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
//	return score(me.score_matrix, (*entry1._seq)[position(entry1)] , (*entry2._seq)[position(entry2)]); // Normal Score using the substitutional matrix
	return score(me.score_matrix, (*entry1._seq)[position(entry1)] ,
                 (*entry2._seq)[position(entry2)]) + me.getMapLineValue(position(entry1), position(entry2));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_ALIGN_RNA_SCORE_STRUCTURE_RNA_H_
