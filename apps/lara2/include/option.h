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
// This file contains the option structures and functions of seqan_laragu
// application.
// ==========================================================================

#ifndef INCLUDE_OPTION_H_
#define INCLUDE_OPTION_H_

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <string.h>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

#include <seqan/arg_parse.h>

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// defines all the constants used in the app
#include "data_types.h"

//using namespace std;

// ============================================================================
// Functors
// ============================================================================

typedef EqualsChar<'.'>        IsDot;
typedef EqualsChar<'/'>        IsSlash;
typedef EqualsChar<'\\'>       IsBackSlash;

#ifdef PLATFORM_WINDOWS
typedef IsBackSlash        IsPathDelimited;
#else
typedef IsSlash            IsPathDelimited;
#endif

// ----------------------------------------------------------------------------
// Class Options
// ----------------------------------------------------------------------------
//TODO this can be re-written as a set of three maps with keys and values but
//if in the future should be embedded in the parsing function is better to
//conserve this style
struct LaraOptions
{
// use base pairs or structure
    bool useBasePairs;
// Name of input file
    seqan::CharString inFile;
// Name of input fileRef
    seqan::CharString inFileRef;
// Name of output file (default: stdout)
    seqan::CharString outFile;
// temporary directory where to save intermediate files. Default: use the input file directory.
    seqan::CharString tmpDir;
// Use the gap scheme to be used in the N-W alignment (default: affine(0))
    unsigned affineLinearDgs;
// Use the global local or global-Unconstrained algorithm (default: global(0) - local(1) )
    bool globalLocal;
// type used for the global-unconstrained alignment AlignConfig <TTop, Tleft, TRight, TDown>
    bool unTop;
    bool unLeft;
    bool unRight;
    bool unDown;
// Threshold of ratio sizes for automatic choice global and global-Unconstrained algorithm (default: 2/3)
    double thrGlobalUnconstr;
// Threshold of ratio sizes for automatic choice global and Local algorithm (default: 1/2))
    double thrGlobalLocal;
// Parameter used during the RNAfold execution to select the minimum energy to be considered
    double thrBppm;
// number of iterations
    unsigned iterations;
// number of non-decreasing iterations
    unsigned nonDecreasingIterations;
// method to be used for the computation of the Lower bound (MWM or approximation can be chosen)
    unsigned lowerBoundMethod;
// value to be considered for the equality of upper and lower bounds difference
    double epsilon;
// my, necessary for computing appropriate step sizes
    double my; //FIXME to be changed the name
// scoring matrix name that should be used for scoring alignment edges in the actual problem
    seqan::CharString laraScoreMatrixName;
//    Score<double, ScoreMatrix<Rna5, Default> > laraScoreMatrix;
    TScoreMatrix laraScoreMatrix;
//    TScoringSchemeRib laraScoreMatrixRib;
// Gap open and extend costs for generating the alignment edges
    double generatorGapOpen;
    double generatorGapExtend;
    //double generatorSuboptimality;
    // FIXME what means this parameter? (Parameter for the generation of alignment edges. The higher the value of 'generatorsuboptimality', the more alignment edges are created.)
// Gap open and extend costs for generating the alignment edges
    double laraGapOpen;
    double laraGapExtend;
// scaling factor for the scores of the alignment edges
    double sequenceScale; // Specifies the contribution of the sequence scores (specified by the larascore matrix) to the overall structural alignment.
// gap penalty for RSA
    double rsaGapPenalty;
// scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM
    unsigned structureScoring;
// define the weight of _half_ an interaction match for fixed structures
    double fixedStructWeight;
// if structureScoring=SCALING then we have to give a scaling factor
    double scalingFactor;
// if true the opposite lines found during the upper bound computation with a max score must be saved and used
    bool useOppositLineUB;
// specify the location of T-COFFEE
    seqan::CharString tcoffeeLocation;
// specify the method to be used to create the T-Coffe library
    unsigned tcoffeLibMode;
// Define the id of the sequence that must be splitted
    unsigned splitSequence;
// window size specifies the length of the sliding window when the local alignment algorithm is used
// (long sequence vs short)
    unsigned windowSize;
// time used for an hard timeout
    unsigned timeLimit;
// verbose(0) no outputs,
// verbose(1) Displays global statistics,
// verbose(2) Displays extensive statistics for each batch of reads,
// verbose(3) Debug output.
    unsigned verbose;
// number of threads forced
    unsigned threads;
// number of threads detected
    unsigned threadsCount;
    LaraOptions() :
            useBasePairs(true),
            affineLinearDgs(0),
            globalLocal(false),
            unTop(false),
            unLeft(false),
            unRight(false),
            unDown(false),
            thrGlobalUnconstr(0.666667),
            thrGlobalLocal(0.5),
            thrBppm(1e-15), // 0.1 is the value used in the old Lara
            iterations(500),
            nonDecreasingIterations(50u),
            lowerBoundMethod(LBLEMONMWM),
            epsilon(0.0001),
            my(1.0),
            laraScoreMatrixName(""), //laraScoreMatrixName("RIBOSUM65"),
            generatorGapOpen(-6.0),
            generatorGapExtend(-2.0),
            // generatorSuboptimality(40),
            laraGapOpen(-6.0),
            laraGapExtend(-2.0),
            sequenceScale(1.0),
            rsaGapPenalty(3.0),
            structureScoring(LOGARITHMIC),
            fixedStructWeight(8.0),
            scalingFactor(1.0),
            useOppositLineUB(true),
            tcoffeeLocation("t_coffee/t_coffee_5.05"),
            tcoffeLibMode(SWITCH),
            splitSequence(1),
            windowSize(100),
            timeLimit(-1),
            verbose(0)
    {
#ifdef _OPENMP
        threadsCount = std::thread::hardware_concurrency(); // omp_get_num_threads()
        threads = threadsCount;
#else
        threadsCount = 1;
#endif
    }
};

// ============================================================================
// Functions
// ============================================================================

using namespace seqan;
// ----------------------------------------------------------------------------
// Function setupArgumentParser()
// ----------------------------------------------------------------------------
template <typename TOption>
void setupArgumentParser(ArgumentParser & parser, TOption const & /* options */)
{
    setAppName(parser, "LaRA-GU");
    setShortDescription(parser, "Lagrangian Relaxation Structural Alignment Algorithm");
    setCategory(parser, "Structural Alignment Algorithm");
    setVersion(parser, "1.0");
    setDate(parser, "2016");
    //setDateAndVersion(parser);
    //setDescription(parser);
    addUsageLine(parser, "./lara <\\fI-i inFile\\fP> \
            [\\fI-w outFile\\fP] [\\fI -parameters\\fP]");
    addOption(parser, ArgParseOption("v", "verbose", "verbose(0) no outputs, verbose(1) Displays global statistics, "
            "verbose(2) Displays extensive statistics for each batch of reads, verbose(3) Debug output.",
            ArgParseArgument::INTEGER, "INT"));

    addSection(parser, "LaRA Alignment Options");
    addOption(parser, ArgParseOption("s", "useBasePairs", "Use structure prediction or fixed structure from extended "
            "input file."));
    addOption(parser, ArgParseOption("g", "affineLinearDgs", "Chose the gap scheme affine(0) linear(1) or dynamic(2) "
            "to be used in the alignment (default: affine(0)). ", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("a", "globalLocal", "Use the local or global algorithm (default: global(false). "));
    addOption(parser, ArgParseOption("ut", "unTop", "type used for the global-unconstrained alignment AlignConfig TTop "
            "(default: top(false). "));
    addOption(parser, ArgParseOption("ul", "unLeft", "type used for the global-unconstrained alignment AlignConfig "
            "TLeft (default: left(false). "));
    addOption(parser, ArgParseOption("ur", "unRight", "type used for the global-unconstrained alignment AlignConfig "
            "TRight (default: right(false). "));
    addOption(parser, ArgParseOption("ud", "unDown", "type used for the global-unconstrained alignment AlignConfig "
            "TDown (default: dow(false). "));
    addOption(parser, ArgParseOption("tgu", "thrGlobalUnconstr", "Threshold of ratio sizes for automatic choice "
            "global and global-Unconstrained algorithm (default: 2/3)", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("tgl", "thrGlobalLocal", "Threshold of ratio sizes for automatic choice global "
            "and Local algorithm (default: 1/2))", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("tb", "thrBppm", "(Parameter used during the RNAfold execution to select the "
            "minimum energy to be considered (default: 1e-15) ", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("iter", "iterations", "number of iterations. ", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("nditer", "nonDecreasingIterations", "number of non-decreasing iterations.",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("lbm", "lowerBoundMethod", "method to be used for the computation of the Lower "
            "bound (MWM or approximation can be chosen)", ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("ep", "epsilon", "value to be considered for the equality of upper and lower "
            "bounds difference", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("my", "my", "necessary for computing appropriate step sizes.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lsm","laraScoreMatrixName", "scoring matrix name that should be used for scoring "
            "alignment edges in the actual problem", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("ggo", "generatorGapOpen",
                                     "Gap open costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("gge", "generatorGapExtend",
                                     "Gap extend costs for generating the alignment edges.",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    // addOption(parser, ArgParseOption("gso", "generatorSuboptimality",
    //                                  "suboptimality costs for generating the alignment edges.",
    //                                  ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lgo", "laraGapOpen",
                                     "Gap open costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("lge", "laraGapExtend",
                                     "Gap extend costs for generating the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("ssc", "sequenceScale",
                                     "scaling factor for the scores of the alignment edges",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("rsag", "rsaGapPenalty", "gap penalty for RSA",
                                     ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("stsc", "structureScoring", "scoring mode, either LOGARITHMIC, SCALE, ORIGINAL, RIBOSUM",
                                     ArgParseArgument::INTEGER, "INT"));
    addOption(parser, ArgParseOption("fsw", "fixedStructWeight", "define the weight of _half_ an interaction match for "
            "fixed structures", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("scal", "scalingFactor", "if structurescoring=SCALING then we have to give a "
            "scaling factor", ArgParseArgument::DOUBLE, "DOUBLE"));
    addOption(parser, ArgParseOption("uolub", "useOppositLineUB", "if true the opposite lines found during the upper "
            "bound computation with a max score must be saved and used", ArgParseArgument::BOOL, "BOOL"));
    addOption(parser, ArgParseOption("spseq", "splitSequence", "Define the id of the sequence that must be splitted.",
                                     ArgParseArgument::INTEGER, "INT")); // TODO fix the meaning of this parameter
    addOption(parser, ArgParseOption("ws", "windowSize", "window size specifies the length of the sliding window.",
                                     ArgParseArgument::INTEGER, "INT")); // TODO fix the meaning of this parameter
    addOption(parser, ArgParseOption("tl", "timeLimit", "some additional option.", ArgParseArgument::INTEGER, "INT"));


    addSection(parser, "Input Options");
    addOption(parser, ArgParseOption("i", "inFile", "Path to the input file", ArgParseArgument::INPUT_FILE, "IN"));
    addOption(parser, ArgParseOption("ir", "inFileRef", "Path to the reference input file",
                                     ArgParseArgument::INPUT_FILE, "IN"));

    addSection(parser, "Output Options");
    addOption(parser, seqan::ArgParseOption( "w", "outFile", "Path to the output file (default: stdout)",
                                             ArgParseArgument::OUTPUT_FILE, "OUT"));
    addOption(parser, ArgParseOption("td", "tmpDir", "Specify a temporary directory where to save intermediate files. \
            Default: use the input file directory.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("tcl", "tcoffeeLocation", "location of T-COFFEE.", ArgParseOption::STRING));
    addOption(parser, ArgParseOption("tcm","tcoffeLibMode", "method used to create the T-Coffe library "
                                             "either     PROPORTIONAL, SWITCH (defoult), ALLINTER, FIXEDINTER", ArgParseArgument::INTEGER, "INT"));

    // Setup performance options.
    addSection(parser, "Performance Options");
    addOption(parser, ArgParseOption("t", "threads", "Specify the number of threads to use.", ArgParseOption::INTEGER));
    setMinValue(parser, "threads", "1");
#ifdef _OPENMP
    setMaxValue(parser, "threads", std::to_string(std::thread::hardware_concurrency() + 1));
#else
    setMaxValue(parser, "threads", "1");
#endif
//    setDefaultValue(parser, "threads", options.threadsCount);
}

// ----------------------------------------------------------------------------
// Function setEnv()
// ----------------------------------------------------------------------------

template <typename TString, typename TValue>
bool setEnv(TString const & key, TValue const & value)
{
#ifdef PLATFORM_WINDOWS
    return !_putenv_s(toCString(key), toCString(value));
#else
    return !setenv(toCString(key), toCString(value), true);
#endif
}

// ----------------------------------------------------------------------------
// Function getCwd()
// ----------------------------------------------------------------------------

template <typename TString>
void getCwd(TString & string)
{
    char cwd[1000];

#ifdef PLATFORM_WINDOWS
    _getcwd(cwd, 1000);
#else
    ignoreUnusedVariableWarning(getcwd(cwd, 1000));
#endif

    assign(string, cwd);
}

// ----------------------------------------------------------------------------
// Function firstOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
firstOf(TString const & string, TFunctor const & f)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = begin(string, Standard());
    skipUntil(it, f);

    return it;
}

// ----------------------------------------------------------------------------
// Function lastOf()
// ----------------------------------------------------------------------------

template <typename TString, typename TFunctor>
inline typename Iterator<TString const, Standard>::Type
lastOf(TString const & string, TFunctor const & f)
{
    typedef ModifiedString<TString const, ModReverse>        TStringRev;
    typedef typename Iterator<TStringRev, Standard>::Type    TIterRev;

    TStringRev revString(string);
    TIterRev revIt = firstOf(revString, f);

    return end(string) - position(revIt, revString);
}

// ----------------------------------------------------------------------------
// Function getPath()
// ----------------------------------------------------------------------------

template <typename TString>
inline typename Prefix<TString const>::Type
getPath(TString const & string)
{
    typedef typename Iterator<TString const, Standard>::Type TIter;

    TIter it = lastOf(string, IsPathDelimited());

    if (it != begin(string, Standard())) --it;

    return prefix(string, it);
}

// ----------------------------------------------------------------------------
// Function parseCmd()
// ----------------------------------------------------------------------------
template <typename TOption>
ArgumentParser::ParseResult parse(TOption & options, ArgumentParser & parser, int argc, char const ** argv)
{
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
        return res;
    getOptionValue(options.verbose, parser, "verbose");
    getOptionValue(options.useBasePairs, parser, "useBasePairs");
    getOptionValue(options.affineLinearDgs, parser, "affineLinearDgs");
    getOptionValue(options.globalLocal, parser, "globalLocal");
    getOptionValue(options.unTop, parser, "unTop");
    getOptionValue(options.unLeft, parser, "unLeft");
    getOptionValue(options.unRight, parser, "unRight");
    getOptionValue(options.unDown, parser, "unDown");
    getOptionValue(options.thrGlobalUnconstr, parser, "thrGlobalUnconstr");
    getOptionValue(options.thrGlobalLocal, parser, "thrGlobalLocal");
    getOptionValue(options.thrBppm, parser, "thrBppm");
    getOptionValue(options.iterations, parser, "iterations");
    getOptionValue(options.nonDecreasingIterations, parser, "nonDecreasingIterations");
    getOptionValue(options.lowerBoundMethod, parser, "lowerBoundMethod");
    getOptionValue(options.epsilon, parser, "epsilon");
    getOptionValue(options.my, parser, "my");
    getOptionValue(options.laraScoreMatrixName, parser, "laraScoreMatrixName");
    getOptionValue(options.generatorGapOpen, parser, "generatorGapOpen");
    getOptionValue(options.generatorGapExtend, parser, "generatorGapExtend");
    // getOptionValue(options.generatorSuboptimality, parser, "generatorSuboptimality");
    getOptionValue(options.laraGapOpen, parser, "laraGapOpen");
    getOptionValue(options.laraGapExtend, parser, "laraGapExtend");
    getOptionValue(options.sequenceScale, parser, "sequenceScale");
    getOptionValue(options.rsaGapPenalty, parser, "rsaGapPenalty");
    getOptionValue(options.structureScoring, parser, "structureScoring");
    getOptionValue(options.fixedStructWeight, parser, "fixedStructWeight");
    getOptionValue(options.scalingFactor, parser, "scalingFactor");
    getOptionValue(options.useOppositLineUB, parser, "useOppositLineUB");
    getOptionValue(options.tcoffeeLocation, parser, "tcoffeeLocation");
    getOptionValue(options.tcoffeLibMode, parser, "tcoffeLibMode");
    getOptionValue(options.splitSequence, parser, "splitSequence");
    getOptionValue(options.windowSize, parser, "windowSize");
    getOptionValue(options.timeLimit, parser, "timeLimit");
    getOptionValue(options.threads, parser, "threads");

    getOptionValue(options.inFile, parser, "inFile");
    if (empty(options.inFile))
        return ArgumentParser::PARSE_ERROR;

    getOptionValue(options.inFileRef, parser, "inFileRef");
    getOptionValue(options.outFile, parser, "outFile");
    if (isSet(parser, "outFile"))
    {
        _V(options, "The specified output file is " << options.outFile);
    }
    else
    {
        CharString tmpDir;
        getOptionValue(tmpDir, parser, "tmpDir");
        if (!isSet(parser, "tmpDir"))
        {
            tmpDir = SEQAN_TEMP_FILENAME();
            // remove "/test_file" suffix
            erase(tmpDir, length(tmpDir) - 10u, length(tmpDir));
        }
        setEnv("TMPDIR", tmpDir);
        options.tmpDir = tmpDir;
        _V(options, "The absolute path where to create the tmpDir is " << tmpDir);
    }
    setScoreMatrix(options);
//    showScoringMatrix(options.laraScoreMatrix);
    return ArgumentParser::PARSE_OK;
}

#endif /* INCLUDE_OPTION_SRTUCT_FUNCT_H_ */
