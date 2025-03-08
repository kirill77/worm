#include "pch.h"
#include "Organism.h"
#include "Cell.h"

void Organism::simulateStep(double dt)
{
    for (uint32_t u = 0; u < m_pCells.size(); ++u)
    {
        m_pCells[u]->update(dt);
    }
}
