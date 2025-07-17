#pragma once

#include <memory>

struct World
{
    std::shared_ptr<class Medium> m_pMedium;
    std::shared_ptr<struct Organism> m_pOrganism;
    double m_fCurTimeSec;

public:
    World(std::shared_ptr<Organism> pOrganism);
    void simulateStep(double dt);
    
    double getCurrentTime() const { return m_fCurTimeSec; }
};

