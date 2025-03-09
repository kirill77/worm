#include "pch.h"
#include "Chromosome.h"
#include "Cell.h"
#include "Medium.h"
#include "Spindle.h"

void Chromosome::update(double dt, Cell& cell, Medium& medium)
{
    // Update based on cell cycle state
    switch (cell.getCellCycleState())
    {
        case CellCycleState::PROPHASE:
            condense(static_cast<float>(dt));
            break;
            
        case CellCycleState::METAPHASE:
            if (!m_isAttached && !m_isSeparated)
            {
                // Try to attach to spindle if not already attached
                if (auto pSpindle = cell.getSpindle())
                {
                    tryAttachToSpindle(*pSpindle);
                }
            }
            break;
            
        case CellCycleState::ANAPHASE:
            if (m_isAttached && !m_isSeparated)
            {
                // Initiate separation
                separate();
            }
            if (m_isAttached && m_isSeparated)
            {
                // Move along spindle after separation
                if (auto pSpindle = cell.getSpindle())
                {
                    moveAlongSpindle(*pSpindle, static_cast<float>(dt));
                }
            }
            break;
            
        case CellCycleState::TELOPHASE:
            decondense(static_cast<float>(dt));
            break;
    }
}

void Chromosome::condense(float dt)
{
    m_condensation = std::min(1.0f, m_condensation + CONDENSATION_RATE * dt);
}

void Chromosome::decondense(float dt)
{
    m_condensation = std::max(0.0f, m_condensation - DECONDENSATION_RATE * dt);
}

void Chromosome::separate()
{
    if (m_isAttached && !m_isSeparated)
    {
        m_isSeparated = true;
        // Initial separation creates a small gap between chromatids
        m_position = m_position + float3(0.0f, SEPARATION_DISTANCE, 0.0f);
    }
}

bool Chromosome::tryAttachToSpindle(const Spindle& spindle)
{
    if (!m_isAttached && spindle.isAssembled())
    {
        // Calculate distance to spindle
        float3 spindlePos = spindle.getPosition();
        float3 toSpindle = spindlePos - m_position;
        float dist = length(toSpindle);
        
        // Attach if close enough
        if (dist < 0.2f)
        {
            m_isAttached = true;
            m_attachmentPoint = m_position;
            return true;
        }
    }
    return false;
}

void Chromosome::moveAlongSpindle(const Spindle& spindle, float dt)
{
    if (m_isAttached && m_isSeparated)
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
            m_position = m_position + direction * (0.5f * dt);
        }
    }
} 