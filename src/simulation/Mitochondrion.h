#pragma once

#include "Organelle.h"

class Mitochondrion : public Organelle
{
private:
    double m_fNumber;  // Number of mitochondria in the population
    
    static constexpr double ATP_PRODUCTION_RATE = 1.0;
    static constexpr double INITIAL_MITOCHONDRIA = 1e5;  // Start with 100,000 mitochondria

public:
    Mitochondrion()
        : m_fNumber(INITIAL_MITOCHONDRIA) {}

    void update(double dt, Cell& cell, Medium& medium) override;
};

