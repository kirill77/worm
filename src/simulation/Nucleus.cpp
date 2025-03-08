#include "pch.h"
#include "Nucleus.h"
#include "Cell.h"
#include "Medium.h"
#include "MRNA.h"
#include <algorithm>

void Nucleus::update(double dt, CellCycleState cellState, 
                    std::function<void(std::shared_ptr<MRNA>)> addMRNA)
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
        
        // Add mRNAs using the callback
        for (const auto& mRNA : mRNAs)
        {
            addMRNA(mRNA);
        }
    }
}
