#include "pch.h"
#include "Cell.h"
#include "Nucleus.h"
#include "Mitochondrion.h"

Cell::Cell(std::shared_ptr<Medium> pMedium)
{
    m_pMedium = pMedium;

    m_pOrganelles.emplace_back(new Nucleus());
    m_pOrganelles.emplace_back(new Mitochondrion());
    // add other organelles as needed
}
void Cell::simulateStep(double dt)
{
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(dt);
    }
    // Additional interactions or global updates...
}
