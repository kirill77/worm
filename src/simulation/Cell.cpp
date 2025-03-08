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
    // Update all organelles
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(dt, m_cellCycleState, *m_pMedium);
    }
    
    // Check for cell cycle transitions based on conditions
    checkCellCycleTransitions();
    
    // Update medium
    m_pMedium->update(dt);
}

void Cell::checkCellCycleTransitions()
{
    // Get key protein concentrations
    float3 center(0, 0, 0);
    double cdk1 = m_pMedium->getProteinNumber("CDK-1", center);
    double cyclinB = m_pMedium->getProteinNumber("CYB-1", center);
    double plk1 = m_pMedium->getProteinNumber("PLK-1", center);
    
    // Check conditions for each transition
    switch (m_cellCycleState)
    {
        case CellCycleState::INTERPHASE:
            // Transition to prophase when CDK1-CyclinB complex reaches threshold
            if (cdk1 > 1000 && cyclinB > 1000)
            {
                m_cellCycleState = CellCycleState::PROPHASE;
            }
            break;

        case CellCycleState::PROPHASE:
            // Transition to metaphase when chromosomes are condensed and aligned
            // TODO: Add chromosome state monitoring
            break;

        case CellCycleState::METAPHASE:
            // Transition to anaphase when spindle checkpoint is satisfied
            // TODO: Add spindle checkpoint monitoring
            break;

        case CellCycleState::ANAPHASE:
            // Transition to telophase when chromosomes reach poles
            // TODO: Add chromosome position monitoring
            break;

        case CellCycleState::TELOPHASE:
            // Transition to cytokinesis when nuclear envelopes reform
            // TODO: Add nuclear envelope monitoring
            break;

        case CellCycleState::CYTOKINESIS:
            // Transition back to interphase when cell division completes
            // TODO: Add cytokinesis completion monitoring
            break;
    }
}
