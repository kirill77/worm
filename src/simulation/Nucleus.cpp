#include "pch.h"
#include "Nucleus.h"
#include "MRNA.h"
#include "Medium.h"
#include "Cell.h"
#include <algorithm>

void Nucleus::update(double dt, CellCycleState cellState, Medium& extMedium)
{
    // 1. Nuclear envelope dynamics
    switch (cellState)
    {
        case CellCycleState::PROPHASE:
            // Nuclear envelope breaks down
            m_envelopeIntegrity = std::max(0.0, m_envelopeIntegrity - dt * ENVELOPE_BREAKDOWN_RATE);
            break;

        case CellCycleState::TELOPHASE:
            // Nuclear envelope reforms
            m_envelopeIntegrity = std::min(1.0, m_envelopeIntegrity + dt * ENVELOPE_REFORM_RATE);
            break;
    }

    // 2. Transcription (only during interphase when envelope is mostly intact)
    if (cellState == CellCycleState::INTERPHASE && m_envelopeIntegrity > 0.8)
    {
        // Transcribe genes
        auto mRNAs = m_pDNA->transcribeAll(dt);
        
        // Add mRNAs to medium near nucleus
        for (const auto& mRNA : mRNAs)
        {
            // Add mRNAs slightly offset from center to simulate nuclear pores
            float angle = static_cast<float>(rand()) / RAND_MAX * 6.28318f;  // Random angle
            float radius = 0.2f;  // Distance from center
            float3 position(
                radius * cos(angle),
                radius * sin(angle),
                0.0f
            );
            extMedium.addMRNA(mRNA, position);
        }
    }
}
