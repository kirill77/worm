#pragma once

#include "Organelle.h"

class Mitochondrion : public Organelle
{
    double m_fAtp; // current ATP energy pool

public:
    void update(double dt) override
    {
        // Simulate ATP production (cellular respiration)
        // e.g., increase ATP based on available resources
    }
};

