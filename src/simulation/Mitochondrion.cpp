#include "pch.h"
#include "Mitochondrion.h"
#include "Cell.h"
#include <algorithm>
#include <random>

Mitochondrion::Mitochondrion()
    : m_fNumber(N_INITIAL_MITOCHONDRIA) 
{
    // Initialize random number generator with a seed
    std::random_device rd;
    m_rng.seed(rd());
}

void Mitochondrion::update(double dt, Cell& cell, Medium& medium)
{
    auto cellState = cell.getCellCycleState();

    // ATP Production (proportional to number of mitochondria)
    double fATPToCreate = dt * ATP_PRODUCTION_RATE * m_fNumber;
    
    static constexpr uint32_t numLocations = 16;
    generateATP(medium, numLocations, fATPToCreate / numLocations);

#if 0
    // Mitochondrial division/fusion based on energy needs
    float3 position(0.0f, 0.0f, 0.0f);  // Center position for checking ATP level
    double localATP = medium.getAvailableATP(position);
    if (localATP < Medium::MAX_ATP_PER_CELL * 0.3)  // If ATP is low, increase mitochondria
    {
        m_fNumber *= (1.0 + 0.1 * dt);  // 10% increase per second
    }
    else if (localATP > Medium::MAX_ATP_PER_CELL * 0.7)  // If ATP is high, decrease mitochondria
    {
        m_fNumber *= (1.0 - 0.05 * dt);  // 5% decrease per second
    }
#endif

    // During cell division, prepare to split mitochondria population
    if (cellState == CellCycleState::CYTOKINESIS)
    {
        // Mitochondria will be split between daughter cells
        // This will be handled by the Cell class during actual division
    }
}

void Mitochondrion::generateATP(Medium& medium, uint32_t n, double amount)
{
    // Add ATP to n random positions
    for (uint32_t i = 0; i < n; ++i)
    {
        float3 position = generateRandomPosition();
        medium.addATP(amount, position);
    }
}

float3 Mitochondrion::generateRandomPosition()
{
    // Generate random position within the normalized range [-1, 1]
    std::uniform_real_distribution<float> dist(-1.f, 1.f);
    
    float x = dist(m_rng);
    float y = dist(m_rng);
    float z = dist(m_rng);
    
    return float3(x, y, z);
}
