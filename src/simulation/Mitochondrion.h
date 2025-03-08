#pragma once

#include "Organelle.h"

class Mitochondrion : public Organelle
{
private:
    double m_fNumber;  // Number of mitochondria in the population
    double m_fAtp;     // Current ATP energy pool
    
    static constexpr double ATP_PRODUCTION_RATE = 1.0;
    static constexpr double ATP_CONSUMPTION_RATE = 0.2;
    static constexpr double MAX_ATP = 100.0;
    static constexpr double INITIAL_MITOCHONDRIA = 1e5;  // Start with 100,000 mitochondria

public:
    Mitochondrion()
        : m_fNumber(INITIAL_MITOCHONDRIA)
        , m_fAtp(MAX_ATP / 2) {}

    void update(double dt, CellCycleState cellState, Medium& medium) override;
};

