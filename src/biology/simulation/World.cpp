#include "World.h"
#include "biology/simulation/Organism.h"

World::World(std::shared_ptr<Organism> pOrganism)
{
    m_pOrganism = pOrganism;
}

void World::simulateStep(double dt)
{
    m_timeContext.m_deltaTSec = dt;
    m_timeContext.m_curTSec += dt;
    m_pOrganism->simulateStep(m_timeContext);
}