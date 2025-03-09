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

    // ATP Production (proportional to number of mitochondria)
    double productionRate = ATP_PRODUCTION_RATE * (m_fNumber / INITIAL_MITOCHONDRIA);
    
    // Add ATP to the medium at the cell's location (assumed to be at center for now)
    float3 position(0.0f, 0.0f, 0.0f);  // Center position
    medium.addATP(dt * productionRate, position);

    // Mitochondrial division/fusion based on energy needs
    double localATP = medium.getAvailableATP(position);
    if (localATP < Medium::MAX_ATP_PER_CELL * 0.3)  // If ATP is low, increase mitochondria
    {
        m_fNumber *= (1.0 + 0.1 * dt);  // 10% increase per second
    }
    else if (localATP > Medium::MAX_ATP_PER_CELL * 0.7)  // If ATP is high, decrease mitochondria
    {
        m_fNumber *= (1.0 - 0.05 * dt);  // 5% decrease per second
    }

    // During cell division, prepare to split mitochondria population
    if (cellState == CellCycleState::CYTOKINESIS)
    {
        // Mitochondria will be split between daughter cells
        // This will be handled by the Cell class during actual division
    }
}
