#include "CellSim.h"
#include "biology/organelles/Cell.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
}

void CellSim::update(double dt)
{
    m_pCell->update(dt);
}
