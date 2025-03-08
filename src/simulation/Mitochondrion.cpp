#include "pch.h"
#include "Mitochondrion.h"
#include <algorithm>

// Constants for ATP metabolism
static constexpr double MAX_ATP = 100.0;
static constexpr double ATP_PRODUCTION_RATE = 1.0;
static constexpr double ATP_CONSUMPTION_RATE = 0.2;

void Mitochondrion::update(double dt, CellCycleState cellState, Medium& extMedium)
{
    // 1. ATP Production
    m_fAtp = std::min(MAX_ATP, m_fAtp + dt * ATP_PRODUCTION_RATE);

    // 2. ATP Consumption based on cell state
    double consumptionMultiplier = 1.0;
    
    switch (cellState)
    {
        case CellCycleState::PROPHASE:
        case CellCycleState::METAPHASE:
            consumptionMultiplier = 2.0;  // Higher energy needs during division
            break;
        case CellCycleState::ANAPHASE:
            consumptionMultiplier = 3.0;  // Peak energy need during chromosome separation
            break;
        default:
            consumptionMultiplier = 1.0;
            break;
    }

    m_fAtp = std::max(0.0, m_fAtp - dt * ATP_CONSUMPTION_RATE * consumptionMultiplier);
}
