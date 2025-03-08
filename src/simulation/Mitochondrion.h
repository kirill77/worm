#pragma once

#include "Organelle.h"
#include "Cell.h"
#include <memory>

class Mitochondrion : public Organelle
{
private:
    double m_fAtp;  // current ATP energy pool
    float3 m_position;  // Position relative to cell center
    
    static constexpr double ATP_PRODUCTION_RATE = 1.0;
    static constexpr double ATP_CONSUMPTION_RATE = 0.2;
    static constexpr double MAX_ATP = 100.0;

public:
    Mitochondrion(const float3& initialPos)
        : m_fAtp(MAX_ATP / 2), m_position(initialPos) {}

    void update(double dt, CellCycleState cellState, Medium& extMedium) override;

private:
    void updatePosition(double dt, CellCycleState state)
    {
        // During division, mitochondria should distribute to both future cells
        if (state == CellCycleState::METAPHASE || state == CellCycleState::ANAPHASE)
        {
            // TODO: Implement mitochondrial segregation
            // This will involve moving mitochondria towards the poles
            // based on spindle position and cell division plane
        }
    }
};

