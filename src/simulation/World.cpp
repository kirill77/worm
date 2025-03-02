#include "pch.h"
#include "World.h"
#include "Worm.h"

World::World()
{
    m_pOrganism = std::make_shared<Worm>();
}
void World::simulateStep(double dt)
{
    m_pOrganism->simulateStep(dt);
}