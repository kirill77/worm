#include "pch.h"
#include "Mitochondrion.h"
#include "Cell.h"
#include <algorithm>

// Constants for ATP metabolism
static constexpr double MAX_ATP = 100.0;
static constexpr double ATP_PRODUCTION_RATE = 1.0;
static constexpr double ATP_CONSUMPTION_RATE = 0.2;

void Mitochondrion::update(double dt, Cell& cell, Medium& medium)
{
    auto cellState = cell.getCellCycleState();

    // 1. ATP Production (proportional to number of mitochondria)
    double productionRate = ATP_PRODUCTION_RATE * (m_fNumber / INITIAL_MITOCHONDRIA);
    m_fAtp = std::min(MAX_ATP, m_fAtp + dt * productionRate);

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

    // 3. Mitochondrial division/fusion based on energy needs
    if (m_fAtp < MAX_ATP * 0.3)  // If ATP is low, increase mitochondria
    {
        m_fNumber *= (1.0 + 0.1 * dt);  // 10% increase per second
    }
    else if (m_fAtp > MAX_ATP * 0.7)  // If ATP is high, decrease mitochondria
    {
        m_fNumber *= (1.0 - 0.05 * dt);  // 5% decrease per second
    }

    // 4. During cell division, prepare to split mitochondria population
    if (cellState == CellCycleState::CYTOKINESIS)
    {
        // Mitochondria will be split between daughter cells
        // This will be handled by the Cell class during actual division
    }
}
