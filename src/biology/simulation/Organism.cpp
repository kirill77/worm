#include "Organism.h"
#include "CellSim.h"

void Organism::simulateStep(double dt)
{
    for (uint32_t u = 0; u < m_pCells.size(); ++u)
    {
        m_pCells[u]->update(dt);
    }
}
