#pragma once

#include <memory>

class World
{
    std::shared_ptr<class Medium> m_pMedium;
    std::shared_ptr<class Organism> m_pOrganism;
    double m_fCurTimeSec;

public:
    World(std::shared_ptr<Organism> pOrganism);
    void simulateStep(double dt);
    
    double getCurrentTime() const { return m_fCurTimeSec; }
};

