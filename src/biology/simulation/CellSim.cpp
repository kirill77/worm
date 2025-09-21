#include "CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"

CellSim::CellSim(std::shared_ptr<Cell> pCell)
    : m_pCell(pCell)
{
}

void CellSim::update(const TimeContext& time)
{
    m_pCell->update(time.m_deltaTSec);
}
