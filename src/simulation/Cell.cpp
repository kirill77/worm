#include "pch.h"
#include "Cell.h"
#include "Nucleus.h"
#include "Mitochondrion.h"
#include "MRNA.h"

Cell::Cell(std::shared_ptr<Medium> pMedium)
    : m_pMedium(pMedium)
    , m_cellCycleState(CellCycleState::INTERPHASE)
{
    // Create DNA for the nucleus
    auto pDNA = std::make_shared<DNA>();
    
    // Create organelles
    m_pOrganelles.push_back(std::make_shared<Nucleus>(pDNA));
    m_pOrganelles.push_back(std::make_shared<Mitochondrion>(float3(0, 0, 0)));
    // add other organelles as needed
}

void Cell::update(double dt)
{
    // Create mRNA callback
    auto addMRNA = [this](std::shared_ptr<MRNA> mRNA) {
        // For now, just add mRNAs to the center of the cell
        m_pMedium->addMRNA(mRNA, float3(0, 0, 0));
    };

    // Update all organelles
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(dt, m_cellCycleState, addMRNA);
    }
    
    // Update medium
    m_pMedium->update(dt);
    
    // Additional cell-level updates can go here
}
