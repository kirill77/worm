#include "pch.h"
#include "Nucleus.h"
#include "molecules/MRNA.h"
#include "Medium.h"
#include "Cell.h"
#include <algorithm>

void Nucleus::update(double fDt, Cell& cell, Medium& medium)
{
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
            break;

        case CellCycleState::TELOPHASE:
            // Nuclear envelope reforms during telophase
            m_fEnvelopeIntegrity = std::min(1.0, m_fEnvelopeIntegrity + fENVELOPE_REFORM_RATE * fDt);
            break;
    }

    // 2. Transcription (only during interphase when envelope is mostly intact)
    if (cell.getCellCycleState() == CellCycleState::INTERPHASE && m_fEnvelopeIntegrity > 0.8)
    {
        // Transcribe genes
        auto mRNAs = transcribeAll(fDt);
        
        // Add mRNAs to medium near nucleus (if we have ATP for synthesis)
        for (const auto& mRNA : mRNAs)
        {
            if (cell.consumeATP(ATPCosts::fMRNA_SYNTHESIS))
            {
                // Add mRNAs slightly offset from center to simulate nuclear pores
                float angle = static_cast<float>(rand()) / RAND_MAX * 6.28318f;  // Random angle
                float radius = 0.2f;  // Distance from center
                float3 position(
                    radius * cos(angle),
                    radius * sin(angle),
                    0.0f
                );
                medium.addMRNA(mRNA, position);
            }
        }
    }
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
        // Collect transcripts from all chromosomes
        for (const auto& chromosome : m_chromosomes)
        {
            auto transcripts = chromosome.transcribe(fDt);
            allTranscripts.insert(allTranscripts.end(), transcripts.begin(), transcripts.end());
        }
    }
    
    return allTranscripts;
}
