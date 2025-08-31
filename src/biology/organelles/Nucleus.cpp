#include "pch.h"
#include "Nucleus.h"
#include "chemistry/MRNA.h"
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
            double cytoplasmicLevel = medium.getMoleculeNumber(Molecule(StringDict::idToString(proteinID), ChemicalType::PROTEIN), nucleusCenter);
            double importAmount = cytoplasmicLevel * importRate * fDt;
            
            if (importAmount > 0.0) {
                importProtein(StringDict::idToString(proteinID), importAmount);
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
        auto newMRNAs = transcribeAll(fDt);
        for (auto& mrna : newMRNAs) {
            MRNA& existingMRNA = m_nuclearCompartment.getOrCreateMRNA(mrna->getName());
            
            // If this is a newly created mRNA (with default values), copy the properties
            if (existingMRNA.getNumber() == 0.0) {
                existingMRNA = *mrna;
            } else {
                // Add to existing mRNA count (accumulate transcribed mRNAs)
                existingMRNA.addNumber(mrna->getNumber());
            }
        }
    }
    
    // 4. Try to export existing mRNAs from nuclear pool (if we have ATP and envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5)
    {
        auto& mrnas = m_nuclearCompartment.getMRNAs();
        auto it = mrnas.begin();
        while (it != mrnas.end()) {
            if (cell.consumeATP(ATPCosts::fMRNA_SYNTHESIS)) {
                auto mrnaPtr = std::make_shared<MRNA>(it->second);
                exportMRNA(mrnaPtr);
                it = mrnas.erase(it); // Remove exported mRNA
            } else {
                ++it; // Keep in nucleus, try again next timestep
            }
        }
    }
    
    // 5. Handle mRNA degradation in nuclear pool
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

std::vector<std::shared_ptr<MRNA>> Nucleus::transcribeAll(double fDt) const
{
    std::vector<std::shared_ptr<MRNA>> allTranscripts;
    
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

void Nucleus::importProtein(const std::string& proteinName, double amount)
{
    // Import protein into nuclear compartment (only if envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5 && amount > 0.0) {
        auto& nuclearProteinPop = m_nuclearCompartment.getOrCreateMolPop(Molecule(proteinName, ChemicalType::PROTEIN));
        nuclearProteinPop.m_fNumber += amount;
    }
}

void Nucleus::exportMRNA(std::shared_ptr<MRNA> mRNA)
{
    // Export mRNA to cytoplasm near nucleus (only if envelope is intact)
    if (m_fEnvelopeIntegrity > 0.5) {
        if (auto pCell = getCell()) {
            // Add mRNAs slightly offset from center to simulate nuclear pores
            float angle = static_cast<float>(rand()) / RAND_MAX * 6.28318f;  // Random angle
            float radius = 0.2f;  // Distance from center
            float3 position(
                radius * cos(angle),
                radius * sin(angle),
                0.0f
            );
            pCell->getInternalMedium().addMRNA(mRNA, position);
        }
    }
}
