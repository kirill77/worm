#include "pch.h"
#include "Chromosome.h"
#include "Cell.h"
#include "Medium.h"
#include "Spindle.h"
#include "MRNA.h"

void Chromosome::update(double fDt, Cell& cell, Medium& medium)
{
    // Update based on cell cycle state
    switch (cell.getCellCycleState())
    {
        case CellCycleState::PROPHASE:
            condense(static_cast<float>(fDt));
            break;
            
        case CellCycleState::METAPHASE:
            if (!m_bIsAttached && !m_bIsSeparated)
            {
                // Try to attach to spindle if not already attached
                if (auto pSpindle = cell.getSpindle())
                {
                    tryAttachToSpindle(*pSpindle);
                }
            }
            break;
            
        case CellCycleState::ANAPHASE:
            if (m_bIsAttached && !m_bIsSeparated)
            {
                // Initiate separation
                separate();
            }
            if (m_bIsAttached && m_bIsSeparated)
            {
                // Move along spindle after separation
                if (auto pSpindle = cell.getSpindle())
                {
                    moveAlongSpindle(*pSpindle, static_cast<float>(fDt));
                }
            }
            break;
            
        case CellCycleState::TELOPHASE:
            decondense(static_cast<float>(fDt));
            break;
    }
}

void Chromosome::condense(float fDt)
{
    m_fCondensation = std::min(1.0f, m_fCondensation + fCONDENSATION_RATE * fDt);
}

void Chromosome::decondense(float fDt)
{
    m_fCondensation = std::max(0.0f, m_fCondensation - fDECONDENSATION_RATE * fDt);
}

void Chromosome::separate()
{
    if (m_bIsAttached && !m_bIsSeparated)
    {
        m_bIsSeparated = true;
        // Initial separation creates a small gap between chromatids
        m_position = m_position + float3(0.0f, fSEPARATION_DISTANCE, 0.0f);
    }
}

bool Chromosome::tryAttachToSpindle(const Spindle& spindle)
{
    if (!m_bIsAttached && spindle.isAssembled())
    {
        // Calculate distance to spindle
        float3 spindlePos = spindle.getPosition();
        float3 toSpindle = spindlePos - m_position;
        float dist = length(toSpindle);
        
        // Attach if close enough
        if (dist < 0.2f)
        {
            m_bIsAttached = true;
            m_attachmentPoint = m_position;
            return true;
        }
    }
    return false;
}

void Chromosome::moveAlongSpindle(const Spindle& spindle, float fDt)
{
    if (m_bIsAttached && m_bIsSeparated)
    {
        // Move towards appropriate spindle pole based on position
        float3 targetPole = m_position.y > 0 ? 
            spindle.getPlusPole() : spindle.getMinusPole();
            
        float3 toTarget = targetPole - m_position;
        float dist = length(toTarget);
        
        if (dist > 0.01f)
        {
            // Move towards target pole
            float3 direction = normalize(toTarget);
            m_position = m_position + direction * (0.5f * fDt);
        }
    }
}

std::vector<std::shared_ptr<MRNA>> Chromosome::transcribe(double fDt) const
{
    // Only transcribe when chromosome is not condensed (during interphase)
    if (m_fCondensation < 0.1f && m_pDNA)
    {
        return m_pDNA->transcribeAll(fDt);
    }
    return std::vector<std::shared_ptr<MRNA>>();
} 
