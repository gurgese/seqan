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
    seqan::Rna5String sequence;
};

struct tcoffeeLib
{
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

void computeTCoffeWeightsProportional(tcoffeePair & tcPair, RnaAlignmentTraits & traits)
{
    tcoffeeW tcW{};
    TBestAlign const & bestAlignment = traits.forMinDiff;

    for (int i = length(bestAlignment.lines) - 1; i >= 0 ; --i)
    {
        PositionPair const & line = bestAlignment.lines[i];
        tcW.ntSeqH = line.first + 1;
        tcW.ntSeqV = line.second + 1;

        //if(bestAlignment.weightLineVect[line.second].weight > 0)
        tcW.weight = TCOFFSET;

        //TODO check if this control is necessary
        if (traits.interactions[line.first].count(line.second) > 0)
            tcW.weight += static_cast<int>(TCMULT * traits.interactions[line.first][line.second].weight);

        tcPair.alignWeights.push_back(tcW);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsSwitch()
// ----------------------------------------------------------------------------

void computeTCoffeWeightsSwitch(tcoffeePair & tcPair, RnaAlignmentTraits & traits)
{
    tcoffeeW tcW{};
    TBestAlign const & bestAlignment = traits.forMinDiff;

    for (int i =  length(bestAlignment.lines) - 1; i >= 0 ; --i)
    {
        PositionPair const & line = bestAlignment.lines[i];
        tcW.ntSeqH = line.first + 1;
        tcW.ntSeqV = line.second + 1;

//        if(bestAlignment.weightLineVect[line.second].weight > 0)
        if (traits.interactions[line.first].count(line.second) > 0)
            tcW.weight = TCMAX;
        else
            tcW.weight = TCOFFSET;

        tcPair.alignWeights.push_back(tcW);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsSeqAlignOnly()
// ----------------------------------------------------------------------------

void computeTCoffeWeightsSeqAlignOnly(tcoffeePair & tcPair, RnaAlignmentTraits & traits)
{
    tcoffeeW tcW{};
    TBestAlign const & bestAlignment = traits.forMinDiff;

    for (int i =  length(bestAlignment.lines) - 1; i >= 0 ; --i)
    {
        tcW.ntSeqH = bestAlignment.lines[i].first + 1;
        tcW.ntSeqV = bestAlignment.lines[i].second + 1;

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
    std::size_t size = numVertices(graph.inter);
    resize(flagVect, size, false);
//    std::cout << graph.inter << std::endl;
//    std::cout << "size = " << size << std::endl;
    for (unsigned i = 0; i < size - 1; ++i)
    {
        String<unsigned> adjVect;
        getVertexAdjacencyVector(adjVect, graph.inter, i);
//        std::cout << length(adjVect) << " ";
        for (unsigned j = 0; j < length(adjVect); ++j)
        {
            flagVect[adjVect[j]] = true;
        }
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsAllInter()
// ----------------------------------------------------------------------------

void computeTCoffeWeightsAllInter(tcoffeePair & tcPair, LaraOptions const & options, RnaAlignmentTraits & traits)
{
    tcoffeeW tcW{};
    if (numVertices(traits.bppGraphH.inter) > 1 && numVertices(traits.bppGraphV.inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, traits.bppGraphH);
        createInteractionFlags(flagVectV, traits.bppGraphV);
        TBestAlign const & bestAlignment = traits.forMinDiff;

        for (int i = length(bestAlignment.lines) - 1; i >= 0; --i)
        {
            tcW.ntSeqH = bestAlignment.lines[i].first + 1;
            tcW.ntSeqV = bestAlignment.lines[i].second + 1;

//            std::cout << "pair " << bestAlignment.lines[i].first + 1 << "/"
//                      << bestAlignment.lines[i].second + 1 << std::endl;
//        if(degree(rna1.bppMatrGraphs[0].inter, bestAlignment.lines[i].first) > 0
// && degree(rna2.bppMatrGraphs[0].inter, bestAlignment.lines[i].second) > 0)
            if (flagVectH[bestAlignment.lines[i].first] && flagVectV[bestAlignment.lines[i].second])

                tcW.weight = TCOFFSET;
            else
                tcW.weight = TCMAX;

// if(rna1.fixedGraphs[0].inter[bestAlignment.lines[i].first] == "."
// || rna2.fixedGraphs[0][bestAlignment.lines[i].second] == ".")
// VIENNA NOTATION https://www.tbi.univie.ac.at/RNA/documentation.html
            tcPair.alignWeights.push_back(tcW);
        }
    }
    else
    {
        _VVV(options, " passed from SEQALIGNONLY mode");
        computeTCoffeWeightsSeqAlignOnly(tcPair, traits);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeightsFixedInter()
// ----------------------------------------------------------------------------

void computeTCoffeWeightsFixedInter(tcoffeePair & tcPair, LaraOptions const & options, RnaRecord const & rna1,
                                    RnaRecord const & rna2, RnaAlignmentTraits & traits)
{
    tcoffeeW tcW{};
    if (!empty(rna1.fixedGraphs) && !empty(rna2.fixedGraphs) &&
        numVertices(rna1.fixedGraphs[0].inter) > 1 && numVertices(rna2.fixedGraphs[0].inter) > 1)
    {
        String<bool> flagVectH, flagVectV;
        createInteractionFlags(flagVectH, rna1.fixedGraphs[0]);
        createInteractionFlags(flagVectV, rna2.fixedGraphs[0]);
        TBestAlign const & bestAlignment = traits.forMinDiff;

        for (int i = length(bestAlignment.lines) - 1; i >= 0; --i)
        {
            tcW.ntSeqH = bestAlignment.lines[i].first + 1;
            tcW.ntSeqV = bestAlignment.lines[i].second + 1;

            if (flagVectH[bestAlignment.lines[i].first] && flagVectV[bestAlignment.lines[i].second])
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
        computeTCoffeWeightsSeqAlignOnly(tcPair, traits);
    }
}

// ----------------------------------------------------------------------------
// Function computeTCoffeWeights()
// ----------------------------------------------------------------------------

void computeTCoffeWeights(TTCoffeeLib & tcLib, LaraOptions const & options, RnaStructContentsPair const & filecontents,
                          RnaAlignmentTraitsVector & traitsVect)
{
    #pragma omp parallel for num_threads(options.threads)
    for (RnaAlignmentTraitsVector::iterator traitsIt = traitsVect.begin(); traitsIt < traitsVect.end(); ++traitsIt)
    {
//        RnaAlignmentTraits & traits = alignmentTraits[i];
        tcoffeePair tcPair;
        tcPair.idSeqH = traitsIt->sequenceIndices.first + 1;
        tcPair.idSeqV = traitsIt->sequenceIndices.second + 1;

        RnaRecord const & rna1 = filecontents.first.records[traitsIt->sequenceIndices.first];
        RnaRecord const & rna2 = empty(filecontents.second.records)
                                 ? filecontents.first.records[traitsIt->sequenceIndices.second]
                                 : filecontents.second.records[traitsIt->sequenceIndices.second];

        if (!empty(filecontents.second.records))
            tcPair.idSeqV += filecontents.first.records.size();

        switch (options.tcoffeLibMode)
        {
            case PROPORTIONAL: computeTCoffeWeightsProportional(tcPair, *traitsIt); break;
            case SWITCH:       computeTCoffeWeightsSwitch(tcPair, *traitsIt); break;
            case ALLINTER:     computeTCoffeWeightsAllInter(tcPair, options, *traitsIt); break;
            case FIXEDINTER:   computeTCoffeWeightsFixedInter(tcPair, options, rna1, rna2, *traitsIt); break;
            default:           std::cout << "Select one of the available modes to compute the T-COFFEE library.\n";
        }

        tcLib.rnaPairs.push_back(tcPair);
    }
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
        tcoffeLibFile << tcLib.rnas.size() << std::endl;
        for (unsigned i = 0; i < tcLib.rnas.size(); ++i)
        {
            tcoffeLibFile << tcLib.rnas[i].name << " " << length(tcLib.rnas[i].sequence) << " ";
            tcoffeLibFile <<  tcLib.rnas[i].sequence << std::endl;
        }
        for (unsigned i = 0; i < tcLib.rnaPairs.size(); ++i)
        {
            tcoffeLibFile << "# " << tcLib.rnaPairs[i].idSeqH << " " << tcLib.rnaPairs[i].idSeqV << std::endl;
            for (tcoffeeW const & tcWeight : tcLib.rnaPairs[i].alignWeights)
                tcoffeLibFile << tcWeight.ntSeqH << " " << tcWeight.ntSeqV << " " << tcWeight.weight << std::endl;
        }
        tcoffeLibFile << "! SEQ_1_TO_N" << std::endl;
        tcoffeLibFile.close();
    }
}

void createTCoffeeLib(LaraOptions const & options, RnaStructContentsPair const & filecontents,
                      RnaAlignmentTraitsVector & alignmentTraits)
{
    SEQAN_ASSERT(!empty(alignmentTraits));

    TTCoffeeLib tcLib{};
    tcLib.rnas.reserve(filecontents.first.records.size() + filecontents.second.records.size());

    for (RnaRecord const & record : filecontents.first.records)
        tcLib.rnas.push_back(rnaSeqs{record.name, record.sequence});

    for (RnaRecord const & record : filecontents.second.records)
        tcLib.rnas.push_back(rnaSeqs{record.name, record.sequence});

    computeTCoffeWeights(tcLib, options, filecontents, alignmentTraits);
    writeFileTCoffeeLib(tcLib, options);
}

#endif //_INCLUDE_TCOFFEE_INTERFACE_H_
