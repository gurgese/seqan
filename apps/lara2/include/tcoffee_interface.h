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
// This file contains functions to interface the output of seqan_laragu
// application with the input of T-Coffee app.
// ==========================================================================

#ifndef _INCLUDE_TCOFFEE_INTERFACE_H_
#define _INCLUDE_TCOFFEE_INTERFACE_H_

// ============================================================================
// Prerequisites
// ============================================================================

using namespace seqan;

struct tcoffeeW
{
    unsigned ntSeqH;
    unsigned ntSeqV;
    double weight;
};

struct tcoffeePair
{
    unsigned idSeqH;
    unsigned idSeqV;
    std::vector<tcoffeeW> alignWeights;
};

struct rnaSeqs
{
    seqan::CharString name;
    unsigned length;
    seqan::Rna5String sequence;
};

struct tcoffeeLib
{
    unsigned size;
    std::vector<rnaSeqs> rnas;
    std::vector<tcoffeePair> rnaPairs;
};

typedef tcoffeeLib TTCoffeeLib;

// ----------------------------------------------------------------------------
// std headers
// ----------------------------------------------------------------------------

#include <iostream>
#include <fstream>

// ----------------------------------------------------------------------------
// SeqAn headers
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// App headers
// ----------------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsProportional()
// ----------------------------------------------------------------------------

template <typename TOption>
void computeTCoffeWeightsProportional(tcoffeePair & tcPair, SEQAN_UNUSED TOption const & options, RnaRecord const & rna1,
                                      RnaRecord const & rna2, TRnaAlign & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for(int i =  rnaAlign.forMinBound.maskIndex - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].second + 1;
        }
        if(rnaAlign.forMinBound.upperBoundVect[rnaAlign.forMinBound.mask[i].second].maxProbScoreLine > 0)
        {
            tcW.weight = TCOFFSET + static_cast<int>(TCMULT * rnaAlign.forMinBound.upperBoundVect[rnaAlign.forMinBound.mask[i].second].maxProbScoreLine);
        }
        else
        {
            tcW.weight = TCOFFSET;
        }
        tcPair.alignWeights.push_back(tcW);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsSwitch()
// ----------------------------------------------------------------------------

template <typename TOption>
void computeTCoffeWeightsSwitch(tcoffeePair & tcPair, SEQAN_UNUSED TOption const & options, RnaRecord const & rna1,
                                      RnaRecord const & rna2, TRnaAlign & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for(int i =  rnaAlign.forMinBound.maskIndex - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].second + 1;
        }
        if(rnaAlign.forMinBound.upperBoundVect[rnaAlign.forMinBound.mask[i].second].maxProbScoreLine > 0)
        {
            tcW.weight = TCMAX;
        }
        else
        {
            tcW.weight = TCOFFSET;
        }
        tcPair.alignWeights.push_back(tcW);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsSeqAlignOnly()
// ----------------------------------------------------------------------------

template <typename TOption>
void computeTCoffeWeightsSeqAlignOnly(tcoffeePair & tcPair, SEQAN_UNUSED TOption const & options, RnaRecord const & rna1,
                                      RnaRecord const & rna2, TRnaAlign & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for(int i =  rnaAlign.forMinBound.maskIndex - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.mask[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.mask[i].second + 1;
        }
        tcW.weight = TCOFFSET;
        tcPair.alignWeights.push_back(tcW);
    }
}

// ----------------------------------------------------------------------------
// Function createInteractionFlags()
// ----------------------------------------------------------------------------

void createInteractionFlags(String<bool> & flagVect, RnaStructureGraph const & graph)
{
    clear(flagVect);
    unsigned size = numVertices(graph.inter);
    resize(flagVect, size, false);
//    std::cout << graph.inter << std::endl;
//    std::cout << "size = " << size << std::endl;
    for(unsigned i = 0; i < size - 1; ++i)
    {
        String<unsigned> adjVect;
        getVertexAdjacencyVector(adjVect, graph.inter, i);
//        std::cout << length(adjVect) << " ";
        for(unsigned j = 0; j < length(adjVect); ++j)
        {
            flagVect[adjVect[j]] = true;
        }
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsAllInter()
// ----------------------------------------------------------------------------

template <typename TOption>
void computeTCoffeWeightsAllInter(tcoffeePair & tcPair, TOption const & options, RnaRecord const & rna1,
                                RnaRecord const & rna2, TRnaAlign & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    if(numVertices(rnaAlign.bppGraphH.inter) > 1 && numVertices(rnaAlign.bppGraphV.inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, rnaAlign.bppGraphH);
        createInteractionFlags(flagVectV, rnaAlign.bppGraphV);
        for (int i = rnaAlign.forMinBound.maskIndex - 1; i >= 0; --i)
        {
            if(length(rna1.sequence) >= length(rna2.sequence))
            {
                tcW.ntSeqH = rnaAlign.forMinBound.mask[i].first + 1;
                tcW.ntSeqV = rnaAlign.forMinBound.mask[i].second + 1;
            }
            else
            {
                tcW.ntSeqV = rnaAlign.forMinBound.mask[i].first + 1;
                tcW.ntSeqH = rnaAlign.forMinBound.mask[i].second + 1;
            }
//            std::cout << "pair " << rnaAlign.forMinBound.mask[i].first + 1 << "/"
//                      << rnaAlign.forMinBound.mask[i].second + 1 << std::endl;
            //        if(degree(rna1.bppMatrGraphs[0].inter, rnaAlign.forMinBound.mask[i].first) > 0 && degree(rna2.bppMatrGraphs[0].inter, rnaAlign.forMinBound.mask[i].second) > 0)
            if (flagVectH[rnaAlign.forMinBound.mask[i].first] && flagVectV[rnaAlign.forMinBound.mask[i].second])

                //        if(rna1.fixedGraphs[0].inter[rnaAlign.forMinBound.mask[i].first] == "." || rna2.fixedGraphs[0][rnaAlign.forMinBound.mask[i].second] == ".") // VIENNA NOTATION https://www.tbi.univie.ac.at/RNA/documentation.html
            {
                tcW.weight = TCOFFSET;
            }
            else
            {
                tcW.weight = TCMAX;
            }
            tcPair.alignWeights.push_back(tcW);
        }
    }
    else
    {
        _VV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, options, rna1, rna2, rnaAlign);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsFixedInter()
// ----------------------------------------------------------------------------

template <typename TOption>
void computeTCoffeWeightsFixedInter(tcoffeePair & tcPair, TOption const & options, RnaRecord const & rna1,
                                  RnaRecord const & rna2, TRnaAlign & rnaAlign)
{
    tcoffeeW tcW;
    if(numVertices(rna1.fixedGraphs[0].inter) > 1 && numVertices(rna2.fixedGraphs[0].inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, rna1.fixedGraphs[0]);
        createInteractionFlags(flagVectV, rna2.fixedGraphs[0]);
        for (int i = rnaAlign.forMinBound.maskIndex - 1; i >= 0; --i)
        {
            if(length(rna1.sequence) >= length(rna2.sequence))
            {
                tcW.ntSeqH = rnaAlign.forMinBound.mask[i].first + 1;
                tcW.ntSeqV = rnaAlign.forMinBound.mask[i].second + 1;
            }
            else
            {
                tcW.ntSeqV = rnaAlign.forMinBound.mask[i].first + 1;
                tcW.ntSeqH = rnaAlign.forMinBound.mask[i].second + 1;
            }
            if (flagVectH[rnaAlign.forMinBound.mask[i].first] && flagVectV[rnaAlign.forMinBound.mask[i].second])
            {
                tcW.weight = TCOFFSET;
            }
            else
            {
                tcW.weight = TCMAX;
            }
            tcPair.alignWeights.push_back(tcW);
        }
    }
    else
    {
        _VV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, options, rna1, rna2, rnaAlign);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsMethodSel()
// ----------------------------------------------------------------------------
template <typename TOption, typename TRecord>
void computeTCoffeWeightsMethodSel(tcoffeePair & tcPair, TOption const & options, TRecord const & rna1,
                         TRecord const & rna2, TRnaAlign & rnaAlign)
{
    if(length(rna1.fixedGraphs[0]) > 0 &&
       length(rna2.fixedGraphs[0]) > 0)
    {
        if (options.tcoffeLibMode == PROPORTIONAL) {
            _VV(options, " passed from PROPORTIONAL mode");
            computeTCoffeWeightsProportional(tcPair, options, rna1, rna2, rnaAlign);
        } else if (options.tcoffeLibMode == SWITCH) {
            _VV(options, " passed from SWITCH mode");
            computeTCoffeWeightsSwitch(tcPair, options, rna1, rna2, rnaAlign);
        } else if (options.tcoffeLibMode == ALLINTER) {
//                std::cout << i << " seq1 " << rnaAligns[i].idBppSeqH << " seq2 " << rnaAligns[i].idBppSeqV << std::endl;
            _VV(options, " passed from ALLINTER mode");
            computeTCoffeWeightsAllInter(tcPair, options, rna1, rna2, rnaAlign);
        } else if (options.tcoffeLibMode == FIXEDINTER) {
//                std::cout << i << " seq1 " << rnaAligns[i].idBppSeqH << " seq2 " << rnaAligns[i].idBppSeqV << std::endl;
            _VV(options, " passed from FIXEDINTER mode");
            computeTCoffeWeightsFixedInter(tcPair, options, rna1, rna2, rnaAlign);
        } else {
            std::cout << "Select one of the available modes to compute the T-COFFE library" << std::endl;
            //            return -1;
        }
    }
    else
    {
        _VV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, options, rna1, rna2, rnaAlign);
    }
}
// ----------------------------------------------------------------------------
// Function computeTCoffeWeights()
// ----------------------------------------------------------------------------

template <typename TOption>
int computeTCoffeWeights(TTCoffeeLib & tcLib, TOption const & options, RnaStructContents const & filecontents1,
                          RnaStructContents const & filecontents2, TRnaAlignVect & rnaAligns)
{
//    #pragma omp parallel for num_threads(options.threads)
    for(unsigned i = 0; i < rnaAligns.size(); ++i)
    {
        tcoffeePair tcPair;
        tcPair.idSeqH = rnaAligns[i].idBppSeqH + 1;
        tcPair.idSeqV = filecontents1.records.size() + rnaAligns[i].idBppSeqV + 1;
//        if(rnaAligns[i].forMinBound.upperBound > 0)
        computeTCoffeWeightsMethodSel(tcPair, options, filecontents1.records[rnaAligns[i].idBppSeqH],
                                          filecontents2.records[rnaAligns[i].idBppSeqV], rnaAligns[i]);
        tcLib.rnaPairs.push_back(tcPair);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeights()
// ----------------------------------------------------------------------------

template <typename TOption>
int computeTCoffeWeights(TTCoffeeLib & tcLib, TOption const & options, RnaStructContents const & filecontents1,
                          TRnaAlignVect & rnaAligns)
{
//    #pragma omp parallel for num_threads(options.threads)
    for(unsigned i = 0; i < rnaAligns.size(); ++i)
    {
        tcoffeePair tcPair;
        tcPair.idSeqH = rnaAligns[i].idBppSeqH + 1;
        tcPair.idSeqV = rnaAligns[i].idBppSeqV + 1;
//        std::cout << tcPair.idSeqH << " " << tcPair.idSeqV << std::endl;
//        if(rnaAligns[i].forMinBound.upperBound > 0)
        computeTCoffeWeightsMethodSel(tcPair, options, filecontents1.records[rnaAligns[i].idBppSeqH],
                       filecontents1.records[rnaAligns[i].idBppSeqV], rnaAligns[i]);

        tcLib.rnaPairs.push_back(tcPair);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function createTCoffeeLib()
// ----------------------------------------------------------------------------

template <typename TOption>
void createTCoffeeLib(TOption const & options, bool const & singleOrDoubleInFile,
                      RnaStructContents const & filecontents1, RnaStructContents const & filecontents2, TRnaAlignVect & rnaAligns)
{
    TTCoffeeLib tcLib;
    tcLib.size = filecontents1.records.size() + filecontents2.records.size();
    tcLib.rnas.resize(tcLib.size);
    unsigned l = 0;
    for(unsigned i = 0; i < filecontents1.records.size(); ++i)
    {
        tcLib.rnas[l].name = filecontents1.records[i].name;
        tcLib.rnas[l].length = length(filecontents1.records[i].sequence);
        tcLib.rnas[l].sequence = filecontents1.records[i].sequence;
//        std::cout << filecontents1.records[i].sequence << std::endl;
        ++l;
    }
    if(!singleOrDoubleInFile)
    {
        computeTCoffeWeights(tcLib, options, filecontents1, rnaAligns);
    }
    else
    {
        for(unsigned i = 0; i < filecontents2.records.size(); ++i)
        {
            tcLib.rnas[l].name = filecontents2.records[i].name;
            tcLib.rnas[l].length = length(filecontents2.records[i].sequence);
            tcLib.rnas[l].sequence = filecontents2.records[i].sequence;
            ++l;
        }
        computeTCoffeWeights(tcLib, options, filecontents1, filecontents2, rnaAligns);
    }
    writeFileTCoffeeLib(tcLib, options);
}


template <typename TOption>
void writeFileTCoffeeLib(TTCoffeeLib & tcLib, TOption const & options)
{
    bool had_err = false;
    std::ofstream tcoffeLibFile;
    if (empty(options.outFile))
    {
        seqan::CharString filePath = options.tmpDir;
        append(filePath, "/tcoffeLara.lib");
        tcoffeLibFile.open(toCString(filePath));
        if (!tcoffeLibFile.is_open())
        {
            std::cerr << "Unable to open the tcoffee lib file for writing: " << filePath << std::endl;
            had_err = true;
        }
    }
    else
    {
        tcoffeLibFile.open(toCString(options.outFile), std::ios::out);
        if (!tcoffeLibFile.is_open())
        {
            std::cerr << "Unable to open the specified output file for writing: " << options.outFile << std::endl;
            had_err = true;
        }
    }

    if (!had_err)
    {
        tcoffeLibFile << "! T-COFFEE_LIB_FORMAT_01" << std::endl;
        tcoffeLibFile << tcLib.size << std::endl;
        for (unsigned i = 0; i < tcLib.size; ++i)
        {
            tcoffeLibFile << tcLib.rnas[i].name << " " << tcLib.rnas[i].length << " ";
            tcoffeLibFile <<  tcLib.rnas[i].sequence << std::endl;
        }
        for (unsigned i = 0; i < tcLib.rnaPairs.size(); ++i)
        {
            tcoffeLibFile << "# " << tcLib.rnaPairs[i].idSeqH << " " << tcLib.rnaPairs[i].idSeqV << std::endl;
            for (unsigned j = 0; j < tcLib.rnaPairs[i].alignWeights.size(); ++j)
            {
                tcoffeLibFile << tcLib.rnaPairs[i].alignWeights[j].ntSeqH << " " << tcLib.rnaPairs[i].alignWeights[j].ntSeqV
                          << " " <<  tcLib.rnaPairs[i].alignWeights[j].weight << std::endl;
            }
        }
        tcoffeLibFile << "! SEQ_1_TO_N" << std::endl;
        tcoffeLibFile.close();
    }
}

#endif //_INCLUDE_TCOFFEE_INTERFACE_H_
