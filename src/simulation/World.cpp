#include "pch.h"
#include "World.h"
#include "Organism.h"

World::World(std::shared_ptr<Organism> pOrganism)
{
    m_pOrganism = pOrganism;
}
void World::simulateStep(double dt)
{
    m_pOrganism->simulateStep(dt);
}