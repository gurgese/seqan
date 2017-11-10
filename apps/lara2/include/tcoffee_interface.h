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

void computeTCoffeWeightsProportional(tcoffeePair & tcPair, RnaRecord const & rna1,
                                      RnaRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for (int i = length(rnaAlign.forMinBound.lines) - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].second + 1;
        }

        //if(rnaAlign.forMinBound.weightLineVect[rnaAlign.forMinBound.lines[i].second].weight > 0)
        if(rnaAlign.interactions[rnaAlign.forMinBound.lines[i].first].count(rnaAlign.forMinBound.lines[i].second) > 0) //TODO check if this control is necessary
        {
            //tcW.weight = TCOFFSET + static_cast<int>(TCMULT * rnaAlign.forMinBound.weightLineVect[rnaAlign.forMinBound.lines[i].second].weight);
            tcW.weight = TCOFFSET + static_cast<int>(TCMULT * rnaAlign.interactions[rnaAlign.forMinBound.lines[i].first][rnaAlign.forMinBound.lines[i].second].weight);
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

void computeTCoffeWeightsSwitch(tcoffeePair & tcPair, RnaRecord const & rna1,
                                      RnaRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for(int i =  length(rnaAlign.forMinBound.lines) - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].second + 1;
        }
//        if(rnaAlign.forMinBound.weightLineVect[rnaAlign.forMinBound.lines[i].second].weight > 0)
        if(rnaAlign.interactions[rnaAlign.forMinBound.lines[i].first].count(rnaAlign.forMinBound.lines[i].second) > 0)
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

void computeTCoffeWeightsSeqAlignOnly(tcoffeePair & tcPair, RnaRecord const & rna1,
                                      RnaRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    for(int i =  length(rnaAlign.forMinBound.lines) - 1; i >= 0 ; --i)
    {
        if(length(rna1.sequence) >= length(rna2.sequence))
        {
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].second + 1;
        }
        else
        {
            tcW.ntSeqV = rnaAlign.forMinBound.lines[i].first + 1;
            tcW.ntSeqH = rnaAlign.forMinBound.lines[i].second + 1;
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

void computeTCoffeWeightsAllInter(tcoffeePair & tcPair, LaraOptions const & options, RnaRecord const & rna1,
                                RnaRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
//    std::cout << rnaAlign.forMinBound.bestAlign;
//    std::cout << rnaAlign.forMinBound.bestAlignScore * options.sequenceScale << std::endl;
    tcoffeeW tcW;
    if(numVertices(rnaAlign.bppGraphH.inter) > 1 && numVertices(rnaAlign.bppGraphV.inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, rnaAlign.bppGraphH);
        createInteractionFlags(flagVectV, rnaAlign.bppGraphV);
        for (int i = length(rnaAlign.forMinBound.lines) - 1; i >= 0; --i)
        {
            if(length(rna1.sequence) >= length(rna2.sequence))
            {
                tcW.ntSeqH = rnaAlign.forMinBound.lines[i].first + 1;
                tcW.ntSeqV = rnaAlign.forMinBound.lines[i].second + 1;
            }
            else
            {
                tcW.ntSeqV = rnaAlign.forMinBound.lines[i].first + 1;
                tcW.ntSeqH = rnaAlign.forMinBound.lines[i].second + 1;
            }
//            std::cout << "pair " << rnaAlign.forMinBound.lines[i].first + 1 << "/"
//                      << rnaAlign.forMinBound.lines[i].second + 1 << std::endl;
            //        if(degree(rna1.bppMatrGraphs[0].inter, rnaAlign.forMinBound.lines[i].first) > 0 && degree(rna2.bppMatrGraphs[0].inter, rnaAlign.forMinBound.lines[i].second) > 0)
            if (flagVectH[rnaAlign.forMinBound.lines[i].first] && flagVectV[rnaAlign.forMinBound.lines[i].second])

                //        if(rna1.fixedGraphs[0].inter[rnaAlign.forMinBound.lines[i].first] == "." || rna2.fixedGraphs[0][rnaAlign.forMinBound.lines[i].second] == ".") // VIENNA NOTATION https://www.tbi.univie.ac.at/RNA/documentation.html
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
        _VVV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, rna1, rna2, rnaAlign);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsFixedInter()
// ----------------------------------------------------------------------------

void computeTCoffeWeightsFixedInter(tcoffeePair & tcPair, LaraOptions const & options, RnaRecord const & rna1,
                                    RnaRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
    tcoffeeW tcW;
    if (!empty(rna1.fixedGraphs) && !empty(rna2.fixedGraphs) &&
        numVertices(rna1.fixedGraphs[0].inter) > 1 && numVertices(rna2.fixedGraphs[0].inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, rna1.fixedGraphs[0]);
        createInteractionFlags(flagVectV, rna2.fixedGraphs[0]);
        for (int i = length(rnaAlign.forMinBound.lines) - 1; i >= 0; --i)
        {
            if(length(rna1.sequence) >= length(rna2.sequence))
            {
                tcW.ntSeqH = rnaAlign.forMinBound.lines[i].first + 1;
                tcW.ntSeqV = rnaAlign.forMinBound.lines[i].second + 1;
            }
            else
            {
                tcW.ntSeqV = rnaAlign.forMinBound.lines[i].first + 1;
                tcW.ntSeqH = rnaAlign.forMinBound.lines[i].second + 1;
            }
            if (flagVectH[rnaAlign.forMinBound.lines[i].first] && flagVectV[rnaAlign.forMinBound.lines[i].second])
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
        _VVV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, rna1, rna2, rnaAlign);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsMethodSel()
// ----------------------------------------------------------------------------
template <typename TRecord>
void computeTCoffeWeightsMethodSel(tcoffeePair & tcPair, LaraOptions const & options, TRecord const & rna1,
                         TRecord const & rna2, RnaAlignmentTraits & rnaAlign)
{
    if (options.tcoffeLibMode == PROPORTIONAL) {
        _VVV(options, " passed from PROPORTIONAL mode");
        computeTCoffeWeightsProportional(tcPair, rna1, rna2, rnaAlign);
    } else if (options.tcoffeLibMode == SWITCH) {
        _VVV(options, " passed from SWITCH mode");
        computeTCoffeWeightsSwitch(tcPair, rna1, rna2, rnaAlign);
    } else if (options.tcoffeLibMode == ALLINTER) {
//                std::cout << i << " seq1 " << rnaAligns[i].idBppSeqH << " seq2 " << rnaAligns[i].idBppSeqV << std::endl;
        _VVV(options, " passed from ALLINTER mode");
        computeTCoffeWeightsAllInter(tcPair, options, rna1, rna2, rnaAlign);
    } else if (options.tcoffeLibMode == FIXEDINTER) {
//                std::cout << i << " seq1 " << rnaAligns[i].idBppSeqH << " seq2 " << rnaAligns[i].idBppSeqV << std::endl;
        _VVV(options, " passed from FIXEDINTER mode");
        computeTCoffeWeightsFixedInter(tcPair, options, rna1, rna2, rnaAlign);
    } else {
        std::cout << "Select one of the available modes to compute the T-COFFE library" << std::endl;
        //            return -1;
    }
}
// ----------------------------------------------------------------------------
// Function computeTCoffeWeights()
// ----------------------------------------------------------------------------

int computeTCoffeWeights(TTCoffeeLib & tcLib, LaraOptions const & options,
                         RnaStructContentsPair const & filecontents,
                         RnaAlignmentTraitsVector & rnaAligns)
{
//    #pragma omp parallel for num_threads(options.threads)
    for(unsigned i = 0; i < rnaAligns.size(); ++i)
    {
        tcoffeePair tcPair;
        tcPair.idSeqH = rnaAligns[i].sequenceIndices.first + 1;
        tcPair.idSeqV = filecontents.first.records.size() + rnaAligns[i].sequenceIndices.second + 1;
//        if(rnaAligns[i].forMinBound.upperBound > 0)
        computeTCoffeWeightsMethodSel(tcPair, options, filecontents.first.records[rnaAligns[i].sequenceIndices.first],
                                      filecontents.second.records[rnaAligns[i].sequenceIndices.second], rnaAligns[i]);
        tcLib.rnaPairs.push_back(tcPair);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeights()
// ----------------------------------------------------------------------------

int computeTCoffeWeights(TTCoffeeLib & tcLib, LaraOptions const & options, RnaStructContents const & filecontents1,
                          RnaAlignmentTraitsVector & rnaAligns)
{
//    #pragma omp parallel for num_threads(options.threads)
    for(unsigned i = 0; i < rnaAligns.size(); ++i)
    {
        tcoffeePair tcPair;
        tcPair.idSeqH = rnaAligns[i].sequenceIndices.first + 1;
        tcPair.idSeqV = rnaAligns[i].sequenceIndices.second + 1;
//        std::cout << tcPair.idSeqH << " " << tcPair.idSeqV << std::endl;
//        if(rnaAligns[i].forMinBound.upperBound > 0)
        computeTCoffeWeightsMethodSel(tcPair, options, filecontents1.records[rnaAligns[i].sequenceIndices.first],
                       filecontents1.records[rnaAligns[i].sequenceIndices.second], rnaAligns[i]);

        tcLib.rnaPairs.push_back(tcPair);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function createTCoffeeLib()
// ----------------------------------------------------------------------------

void writeFileTCoffeeLib(TTCoffeeLib & tcLib, LaraOptions const & options)
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

void createTCoffeeLib(LaraOptions const & options, RnaStructContentsPair const & filecontents,
                      RnaAlignmentTraitsVector & rnaAligns)
{
    TTCoffeeLib tcLib;
    tcLib.size = filecontents.first.records.size() + filecontents.second.records.size();
    tcLib.rnas.resize(tcLib.size);
    unsigned l = 0;
    for (unsigned i = 0; i < filecontents.first.records.size(); ++i)
    {
        tcLib.rnas[l].name = filecontents.first.records[i].name;
        tcLib.rnas[l].length = length(filecontents.first.records[i].sequence);
        tcLib.rnas[l].sequence = filecontents.first.records[i].sequence;
//        std::cout << filecontents.first.records[i].sequence << std::endl;
        ++l;
    }

    if (empty(filecontents.second.records))
    {
        computeTCoffeWeights(tcLib, options, filecontents.first, rnaAligns);
    }
    else
    {
        for (unsigned i = 0; i < filecontents.second.records.size(); ++i)
        {
            tcLib.rnas[l].name = filecontents.second.records[i].name;
            tcLib.rnas[l].length = length(filecontents.second.records[i].sequence);
            tcLib.rnas[l].sequence = filecontents.second.records[i].sequence;
            ++l;
        }
        computeTCoffeWeights(tcLib, options, filecontents, rnaAligns);
    }
    writeFileTCoffeeLib(tcLib, options);
}



#endif //_INCLUDE_TCOFFEE_INTERFACE_H_
