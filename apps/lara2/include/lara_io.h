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
// This file contains functions to read the input and to write the output
// of seqan_laragu application.
// ==========================================================================

#ifndef _INCLUDE_LARA_IO_H_
#define _INCLUDE_LARA_IO_H_

// ============================================================================
// Prerequisites
// ============================================================================

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <limits>

#include <seqan/rna_io.h>

using namespace seqan;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRnaFile()
// ----------------------------------------------------------------------------

void readRnaFile(RnaStructContents & filecontents, CharString filename, LaraOptions const & options)
{
    if (empty(filename))
        return;

    RnaStructFileIn rnaStructFile;
    if (open(rnaStructFile, toCString(filename), OPEN_RDONLY))
    {
        _VVV(options, "Input file is RnaStruct.");
        readRecords(filecontents, rnaStructFile, std::numeric_limits<unsigned>::max());
        close(rnaStructFile);
    }
    else
    {
        // Read the file.
        _VVV(options, "Input file is Fasta/Fastq.");
        SeqFileIn seqFileIn(toCString(filename));
        StringSet<CharString> ids;
        StringSet<IupacString> seqs;
        StringSet<CharString> quals;
        readRecords(ids, seqs, quals, seqFileIn);
        close(seqFileIn);

        // Fill the data structures: identifier and sequence.
        resize(filecontents.records, length(ids));
        SEQAN_ASSERT_EQ(length(ids), length(seqs));
        for (typename Size<StringSet<CharString> >::Type idx = 0u; idx < length(ids); ++idx)
        {
            filecontents.records[idx].name = ids[idx];
            filecontents.records[idx].sequence = convert<Rna5String>(seqs[idx]);
        }
        // For FastQ files: add quality annotation.
        if (length(quals) == length(ids))
        {
            for (typename Size<StringSet<CharString> >::Type idx = 0u; idx < length(ids); ++idx)
                filecontents.records[idx].quality = quals[idx];
        }
    }
}

// ----------------------------------------------------------------------------
// Function plotOutput()
// ----------------------------------------------------------------------------

void plotAlignments(LaraOptions const & options, RnaAlignmentTraitsVector & alignmentTraits)
{
    for (RnaAlignmentTraits const & traits : alignmentTraits)
    {
/*
        _VV(options, "******* For Minimum Stepsize *******" << std::endl);
        _VV(options, "Alignment of sequences " << traits.idBppSeqH << ":" << traits.idBppSeqV);
        _VV(options, "Iteration where MinStepSize has been found " << traits.forMinBound.it);
        _VV(options, "Best Lower bound is " << traits.forMinBound.lowerBound);
        _VV(options, "Best Upper bound is " << traits.forMinBound.upperBound);
        _VV(options, "Minimum step size is " << traits.forMinBound.stepSizeBound);
        _VV(options, "Best alignment based on the Min Bounds is " << traits.forMinBound.bestAlignScore << "\n"
                                                                  <<  traits.forMinBound.bestAlign);

        _VV(options, "The step size to be used for Lambda at last iteration is " << traits.stepSize << "\n");

        _VV(options, "+++++++ For Maximum Alignment Score +++++++" << std::endl);
        _VV(options, "Alignment of sequences " << traits.idBppSeqH << ":" << traits.idBppSeqV);
        _VV(options, "Iteration where best score has been found " << traits.forScore.it);
        _VV(options, "Best Lower bound is " << traits.forScore.lowerBound);
        _VV(options, "Best Upper bound is " << traits.forScore.upperBound);
        _VV(options, "Minimum step size is " << traits.forScore.stepSizeBound);
        _VV(options, "The step size to be used for Lambda at last iteration is " << traits.stepSize << "\n");
        _VV(options, "Best alignment based on the Score is " << traits.forScore.bestAlignScore << "\n"
                                                             << traits.forScore.bestAlign);
 */
        _VV(options, "******* For Closest Bounds *******" << std::endl);
        _VV(options, "Alignment of sequences " << traits.sequenceIndices.first << ":" << traits.sequenceIndices.second);
        _VV(options, "Iteration where the closest bounds have been found " << traits.forMinDiff.it);
        _VV(options, "Best Lower bound is " << traits.forMinDiff.lowerBound);
        _VV(options, "Best Upper bound is " << traits.forMinDiff.upperBound);
        _VV(options, "Minimum step size is " << traits.forMinDiff.stepSizeBound);
        _VV(options, "Best alignment based on the the closest Bounds is " << traits.forMinDiff.bestAlignScore << "\n"
                                                                          << traits.forMinDiff.bestAlign);
    }
}

void createFastaAlignmentFile(LaraOptions const & options, RnaAlignmentTraitsVector & alignmentTraits)
{
    std::ofstream fastaFile;
    fastaFile.open(toCString(options.outFile), std::ios::out);
    if (!fastaFile.is_open())
    {
        std::cerr << "Unable to open the specified output file for writing: " << options.outFile << std::endl;
        return;
    }

    RnaAlignment const & ali = alignmentTraits[0].forMinDiff.bestAlign;
    fastaFile << ">sequence1" << std::endl << row(ali, 0) << std::endl;
    fastaFile << ">sequence2" << std::endl << row(ali, 1) << std::endl;
    fastaFile.close();
}

CharString getEbpseqString(RnaStructContents & filecontents)
{
    String<char> outstr;
    if (length(filecontents.records) == 0)
        return outstr;

    RnaIOContext context;
    createPseudoHeader(filecontents.header, filecontents.records);
    writeHeader(outstr, filecontents.header, context, Ebpseq());
    for (RnaRecord & rec : filecontents.records)
        writeRecord(outstr, rec, context, Ebpseq());
    return outstr;
}

#endif //_INCLUDE_LARA_IO_H_
