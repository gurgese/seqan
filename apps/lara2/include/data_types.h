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
// This file contains the variable definition and structures of laragu
// application.
// ==========================================================================

#ifndef _INCLUDE_TOP_DATA_STRUCT_H_
#define _INCLUDE_TOP_DATA_STRUCT_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <utility>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/basic.h>
#include <seqan/rna_io.h>
#include <seqan/align_rna.h>
#include "score_structure_rna.h"

// ============================================================================
// Prerequisites
// ============================================================================

enum LaraFileType
{
    FASTA,
    FASTQ,
    RNASTRUCT,
    UNKNOWN
};

enum LaraScore
{
    LOGARITHMIC,
    SCALE,
    ORIGINAL,
    RIBOSUM
};

enum LaraTCoffeeLibMode
{
    PROPORTIONAL,
    SWITCH,
    ALLINTER,
    FIXEDINTER
};

enum LaraMwmMethod
{
    LBLEMONMWM,
    LBAPPROXMWM,
    LBMWMTEST,
    LBLINEARTIMEMWM
};

// Value to be used to copare the difference between the upper and the lower bound difference
double const EPSILON = 0.0001;

//Values used for T-Coffee lib preparation
int const TCOFFSET = 500;
int const TCMULT = 10000;
int const TCMAX = 1000;

// ============================================================================
// Macro utility
// ============================================================================

#define _V(_opt, _str) {if(_opt.verbose>0) std::cerr << _str << std::endl;}
#define _VV(_opt, _str) {if(_opt.verbose>1) std::cerr << _str << std::endl;}
#define _VVV(_opt, _str) {if(_opt.verbose>2) std::cerr << _str << std::endl;}

// ============================================================================
// Types used in the program
// ============================================================================

typedef seqan::Rna5String TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;      // align type
typedef unsigned TPosition;
typedef double TScoreValue;
typedef seqan::CharString TString;
typedef seqan::Score<double, seqan::ScoreMatrix<seqan::Rna5, seqan::Default> > TScoreMatrix;
typedef seqan::Ribosum65N TRibosum;
typedef float TBioval;
typedef std::map<TPosition, TScoreValue> TMap;
typedef seqan::String<TMap > TMapLine;
typedef std::vector<TMap > TMapVect;
typedef std::vector<seqan::RnaRecord> TRnaVect;
typedef StringSet<Rna5String, Dependent<Generous> > RnaSeqSet;

struct boundStruct
{
    // String with size seq2
    unsigned seq1Index;
    TScoreValue maxProbScoreLine;
    unsigned seq1IndexPairLine;
    unsigned seq2IndexPairLine;
};

// String with size seq2
typedef seqan::String<boundStruct > TBound;

typedef seqan::Graph<seqan::Undirected<double> > TLowerBoundGraph;
//TODO if the Lemon library is used after the tests this graph structure should be chosen as lemon graph in order to
// avoid the copy of the graph

struct lowerBoundLemonStruct
{
    double mwmPrimal; //value of the maximum weighted matching is here saved
    double mwmDual; //value of the maximum weighted matching is here saved
    unsigned mwmCardinality;
};

// String with size seq2
typedef lowerBoundLemonStruct TlowerLemonBound;

// lambda value for subgradient optimization, initialized with 0
struct lambWeightStruct
{
    TScoreValue step;               // actual value of lambda
    TScoreValue maxProbScoreLine;
    unsigned seq1IndexPairLine;
    unsigned seq2IndexPairLine;
};

//typedef std::map<TPosition, lambWeightStruct> TMapWeight;

struct lambStruct
{
    std::map<TPosition, lambWeightStruct> map; //mapLine;
};

//! Lambda Vector
typedef seqan::String<lambStruct> TLambVect;

typedef seqan::Score<double, seqan::ScoreMatrix<seqan::Rna5, TRibosum> > TScoreMatrixRib;
typedef seqan::Score<TScoreValue, RnaStructureScore<TScoreMatrix, TLambVect> > TScoringSchemeStruct;

struct bestAlign
{
    TAlign bestAlign;
    TScoreValue bestAlignScore{std::numeric_limits<TScoreValue>::lowest()};
    int it; //to be used for the best lower bound
    double lowerBound{};
    double upperBound{};
    double stepSizeBound{std::numeric_limits<TScoreValue>::max()};
    TBound upperBoundVect;
    seqan::String<std::pair <unsigned, unsigned> > mask;
    unsigned maskIndex;
};
typedef bestAlign TBestAlign;

struct RnaStructAlign
{
//public:
    seqan::RnaStructureGraph bppGraphH;
    seqan::RnaStructureGraph bppGraphV;
    unsigned idBppSeqH{};
    unsigned idBppSeqV{};
    double sequenceScale{1.0};
// The best computed alignment is saved in these fields
    TBestAlign forScore;
//    TAlign bestAlign;
//    TScoreValue bestAlignScore{std::numeric_limits<TScoreValue>::lowest()};
// Mask that represents the matches from the computed alignment
    seqan::String<unsigned > maskLong;
    seqan::String<std::pair <unsigned, unsigned> > mask;
    unsigned maskIndex;

// Lower bound fields
    double lowerBound{};
//    TBound lowerBoundVect;
// This field is used to approximate the maximum weighted match If tests of this usage are positive we can cosider
// to do not use anymore the Lemon MWM
    TLowerBoundGraph lowerBoundGraph; //graph useful for the seqan::MaximumWeightedMatch() function
    TlowerLemonBound lowerLemonBound;
    double lowerGreedyBound;

// Upper bound fields
    double upperBound{};
    TBound upperBoundVect; // receives interactions of MWM

// Parameters used to compute the stepsize
    int slm{};
    double stepSize{std::numeric_limits<TScoreValue>::max()};
    unsigned nonDecreasingIterations{};
    double my{1.0};

//  Status when the minumum difference between the two bounds is detected
    TBestAlign forMinBound;
//    unsigned itMinBounds; //to be used for the best lower bound
//    double lowerMinBound{};
//    double upperMinBound{};
//    double stepSizeMinBound{std::numeric_limits<TScoreValue>::max()};
//    TAlign bestAlignMinBounds;
//    TScoreValue bestAlignScoreMinBounds;

// String with size seq1 storing all the aligned lines
    TLambVect lamb;

// Scoring scheme used for the structural alignment
    TScoringSchemeStruct structScore;
};// rnaStructAlign;

typedef RnaStructAlign TRnaAlign;
typedef std::vector<TRnaAlign> TRnaAlignVect;

void printRnaStructAlign(TRnaAlign & a, unsigned i)
{
    std::cerr << i << "\tLambda Vector\n";
    unsigned c1 = 0, c2 = 0;
    for (lambStruct & m : a.lamb)
    {
        std::cerr << row(a.forMinBound.bestAlign, 0)[c1++] << " ";
        std::cerr << row(a.forMinBound.bestAlign, 1)[c2++] << " ";
        for (auto p = m.map.begin(); p != m.map.end(); ++p)
        {
            std::cerr << "("<< p->first << "|" << p->second.step << "|" << p->second.maxProbScoreLine << "|" << p->second.seq1IndexPairLine << "|" << p->second.seq2IndexPairLine << ")";
        }
        std::cerr << ";\n";
    }
}

#endif //_INCLUDE_TOP_DATA_STRUCT_H_
