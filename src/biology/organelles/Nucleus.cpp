#include "pch.h"
#include "Nucleus.h"

#include "Medium.h"
#include "Cell.h"
#include <algorithm>

void Nucleus::update(double fDt, Cell& cell)
{
    Medium& medium = cell.getInternalMedium();
    
    // 1. Nuclear import - transport transcription factors from cytoplasm to nucleus
    if (m_fEnvelopeIntegrity > 0.5) {
        float3 nucleusCenter(0.0f, 0.0f, 0.0f);
        
        // List of proteins that should be imported into nucleus (transcription factors)
        static const std::vector<StringDict::ID> transcriptionFactors = {
            StringDict::ID::CDK_2,    // Cell cycle kinase
            StringDict::ID::CCE_1,    // CyclinE
            // Add more transcription factors here as needed
        };
        
        // Simple import model: fraction of cytoplasmic amount enters nucleus per time step
        double importRate = 0.1;  // 10% per time step
        
        for (const auto& proteinID : transcriptionFactors) {
            double cytoplasmicLevel = medium.getMoleculeConcentration(Molecule(proteinID, ChemicalType::PROTEIN), nucleusCenter);
            double importAmount = cytoplasmicLevel * importRate * fDt;

            if (importAmount > 0.0) {
                importMolecule(Molecule(proteinID, ChemicalType::PROTEIN), importAmount);
            }
        }
    }
    
    // Update all chromosomes
    for (auto& chromosome : m_chromosomes)
    {
        chromosome.update(fDt, cell, medium);
    }

    // Update nuclear envelope based on cell cycle state
    switch (cell.getCellCycleState())
    {
        case CellCycleState::PROPHASE:
            // Nuclear envelope breaks down during prophase
            m_fEnvelopeIntegrity = std::max(0.0, m_fEnvelopeIntegrity - fENVELOPE_BREAKDOWN_RATE * fDt);
            // When envelope breaks down, nuclear content mixes with cytoplasm
            if (m_fEnvelopeIntegrity < 0.1) {
                // TODO: Export nuclear contents to cytoplasm during breakdown
            }
            break;

        case CellCycleState::TELOPHASE:
            // Nuclear envelope reforms during telophase
            m_fEnvelopeIntegrity = std::min(1.0, m_fEnvelopeIntegrity + fENVELOPE_REFORM_RATE * fDt);
            break;
    }

    // 3. Transcription and mRNA export (only during interphase when envelope is mostly intact)
    if (cell.getCellCycleState() == CellCycleState::INTERPHASE && m_fEnvelopeIntegrity > 0.8)
    {
        // Transcribe genes using nuclear compartment and add to nuclear pool
        auto newRNAs = transcribeAll(fDt);
        for (auto& rna : newRNAs) {
            Population& existingRNA = m_nuclearCompartment.getOrCreateMolPop(rna->m_molecule);
            
            // Add to existing RNA population (accumulate transcribed RNAs)
            existingRNA.m_fNumber += rna->m_population.m_fNumber;
        }
    }
    
    // 4. Try to export existing RNAs from nuclear pool (if we have ATP and envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5)
    {
        // Find RNA molecules in the nuclear compartment
        auto& molecules = m_nuclearCompartment.m_molecules;
        auto it = molecules.begin();
        while (it != molecules.end()) {
            if (it->first.getType() == ChemicalType::MRNA && it->second.m_fNumber > 0.1) {
                if (cell.consumeATP(ATPCosts::fMRNA_EXPORT)) {
                    auto rnaPtr = std::make_shared<MPopulation>(it->first, it->second.m_fNumber);
                    exportRNA(rnaPtr);
                    it = molecules.erase(it); // Remove exported RNA
                } else {
                    ++it; // Keep in nucleus, try again next timestep
                }
            } else {
                ++it;
            }
        }
    }
    
    // 5. Handle RNA degradation in nuclear pool
    m_nuclearCompartment.updateMRNAs(fDt);
}

bool Nucleus::areChromosomesCondensed() const
{
    for (const auto& chromosome : m_chromosomes)
    {
        if (!chromosome.isFullyCondensed())
        {
            return false;
        }
    }
    return true;
}

bool Nucleus::areChromosomesAttached() const
{
    for (const auto& chromosome : m_chromosomes)
    {
        if (!chromosome.isAttached())
        {
            return false;
        }
    }
    return true;
}

bool Nucleus::areChromosomesSeparated() const
{
    for (const auto& chromosome : m_chromosomes)
    {
        if (!chromosome.isSeparated())
        {
            return false;
        }
    }
    return true;
}

bool Nucleus::areChromosomesDecondensed() const
{
    for (const auto& chromosome : m_chromosomes)
    {
        if (!chromosome.isFullyDecondensed())
        {
            return false;
        }
    }
    return true;
}

std::vector<std::shared_ptr<MPopulation>> Nucleus::transcribeAll(double fDt) const
{
    std::vector<std::shared_ptr<MPopulation>> allTranscripts;
    
    // Only transcribe if nuclear envelope is mostly intact
    if (m_fEnvelopeIntegrity > 0.8)
    {
        // Collect transcripts from all chromosomes using nuclear compartment
        for (const auto& chromosome : m_chromosomes)
        {
            auto transcripts = chromosome.transcribe(fDt, m_nuclearCompartment);
            allTranscripts.insert(allTranscripts.end(), transcripts.begin(), transcripts.end());
        }
    }
    
    return allTranscripts;
}

void Nucleus::importMolecule(const Molecule& molecule, double amount)
{
    // Import molecule into nuclear compartment (only if envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5 && amount > 0.0) {
        auto& nuclearMoleculePop = m_nuclearCompartment.getOrCreateMolPop(molecule);
        nuclearMoleculePop.m_fNumber += amount;
    }
}

void Nucleus::exportRNA(std::shared_ptr<MPopulation> rna)
{
    // Export RNA to cytoplasm near nucleus (only if envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5) {
        if (auto pCell = getCell()) {
            // Add RNAs slightly offset from center to simulate nuclear pores
            float angle = static_cast<float>(rand()) / RAND_MAX * 6.28318f;  // Random angle
            float radius = 0.2f;  // Distance from center
            float3 position(
                radius * cos(angle),
                radius * sin(angle),
                0.0f
            );
            pCell->getInternalMedium().addMolecule(*rna, position);
        }
    }
}
