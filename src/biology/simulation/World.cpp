#include "World.h"
#include "biology/simulation/Organism.h"

World::World(std::shared_ptr<Organism> pOrganism)
    : m_fCurTimeSec(0.0)
{
    m_pOrganism = pOrganism;
}

void World::simulateStep(double dt)
{
    m_fCurTimeSec += dt;
    m_pOrganism->simulateStep(dt);
}