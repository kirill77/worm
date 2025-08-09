#include "Organism.h"
#include "CellSim.h"

void Organism::simulateStep(const TimeContext& time)
{
    for (uint32_t u = 0; u < m_pCellSims.size(); ++u)
    {
        m_pCellSims[u]->update(time);
    }
}
