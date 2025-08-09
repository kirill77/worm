#pragma once

#include <memory>
#include "biology/simulation/TimeContext.h"

class World
{
    std::shared_ptr<class Medium> m_pMedium;
    std::shared_ptr<class Organism> m_pOrganism;
    TimeContext m_timeContext;

public:
    World(std::shared_ptr<Organism> pOrganism);
    void simulateStep(double dt);
    
    double getCurrentTime() const { return m_timeContext.m_curTSec; }
    const TimeContext& getTimeContext() const { return m_timeContext; }
};

