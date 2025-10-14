#include "CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "biology/simulation/PhysicsCore.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
    // Create and initialize physics core; it pulls mesh and volume from the cell
    m_pPhysicsCore = std::make_shared<PhysicsCore>();
    m_pPhysicsCore->initialize(m_pCell);
}

void CellSim::update(const TimeContext& time)
{
    // Advance physics; it pulls volume from cell and pushes mesh back to cortex
    assert(m_pPhysicsCore && m_pCell);
    m_pPhysicsCore->makeTimeStep(time.m_deltaTSec);

    m_pCell->update(time.m_deltaTSec);
}
