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
    MWM_LEMON,
    MWM_GREEDY,
    MWM_SIMPLE

};

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

typedef seqan::Align<seqan::Rna5String, seqan::ArrayGaps> RnaAlignment;      // align type
typedef seqan::Score<double, seqan::ScoreMatrix<seqan::Rna5, seqan::Default> > RnaScoreMatrix;
typedef std::map<unsigned, double> ScMap;
typedef std::vector<std::map<unsigned, double> > InteractionScoreMap;
typedef std::vector<seqan::RnaRecord> RnaRecordVector;
typedef StringSet<Rna5String, Dependent<Generous> > RnaSeqSet;
typedef std::pair<RnaStructContents, RnaStructContents> RnaStructContentsPair;
typedef std::pair<unsigned, unsigned> PositionPair;
typedef seqan::Graph<seqan::Undirected<double> > RnaInteractionGraph;

struct LaraOptions;

// lambda value for subgradient optimization, initialized with 0
struct RnaInteraction
{
    double lambdaValue{0};
    double weight{0};
    double maxProbScoreLine1{0};
    double maxProbScoreLine2{0};
    PositionPair lineL{};
    PositionPair lineM{};
    int iterUpdate{-1};
    bool closedLoop{false};
    // This flag is used for saving mates found during the upper bound update of the line weights.
    bool fromUBPairing{false};
};

typedef seqan::String<std::map<unsigned, RnaInteraction> > OutgoingInteractions;

typedef seqan::Score<double, RnaStructureScore<RnaScoreMatrix, OutgoingInteractions> > LaraScoringScheme;

struct bestAlign
{
    RnaAlignment bestAlign;
    double bestAlignScore{std::numeric_limits<double>::lowest()};
    int it; //to be used for the best lower bound
    double lowerBound{std::numeric_limits<double>::lowest()};
    double upperBound{std::numeric_limits<double>::max()};
    double stepSizeBound{std::numeric_limits<double>::max()};
    //TWeightLine weightLineVect;
    seqan::String<PositionPair> lines;
};
typedef bestAlign TBestAlign;

struct RnaAlignmentTraits
{
    seqan::RnaStructureGraph bppGraphH;
    seqan::RnaStructureGraph bppGraphV;
    std::pair<unsigned, unsigned> sequenceIndices{};

    // Mask that represents the matches from the computed alignment
    seqan::String<PositionPair> lines;

//    TWeightLine lowerBoundVect;
//    This field is used to approximate the maximum weighted match If tests of this usage are positive we can consider
//    to do not use anymore the Lemon MWM
//    TLowerBoundGraph lowerBoundGraph; //graph useful for the seqan::MaximumWeightedMatch() function
//    TlowerLemonBound lowerLemonBound{};

    // Lower bound values (primal)
    double lowerBound{std::numeric_limits<double>::lowest()};
    double bestLowerBound{std::numeric_limits<double>::lowest()};

    // Upper bound values (dual)
    double upperBound{std::numeric_limits<double>::max()};
    double bestUpperBound{std::numeric_limits<double>::max()};

    // Number of edges that violate the relaxed constraint (s_lm).
    int numberOfSubgradients{};
    // Step size for lambda changes (gamma).
    double stepSize{std::numeric_limits<double>::max()};
    // Scaling factor (mu).
    double stepSizeScaling{1.0};
    // Number of iterations without decreasing upper bound.
    unsigned nonDecreasingIterations{0u};
    // Score only from sequence comparison and gap costs.
    double sequenceScore{};

    //  Status when the minimum difference between the two bounds is detected
    TBestAlign forMinBound;
    TBestAlign forMinDiff;
//    TBestAlign forScore;

    // String with size seq1 storing all the aligned lines
    OutgoingInteractions interactions;

    // Scoring scheme used for the structural alignment
    LaraScoringScheme structureScore;
};
typedef std::vector<RnaAlignmentTraits> RnaAlignmentTraitsVector;

#endif //_INCLUDE_TOP_DATA_STRUCT_H_
