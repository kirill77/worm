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
    m_pOrganelles.push_back(std::make_shared<Mitochondrion>());
    // add other organelles as needed
}

void Cell::update(double dt)
{
    // Update all organelles
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(dt, *this, *m_pMedium);
    }
    
    // Check for cell cycle transitions based on conditions
    checkCellCycleTransitions();
    
    // Update medium
    m_pMedium->update(dt);
}

std::shared_ptr<Mitochondrion> Cell::getMitochondrion() const
{
    for (const auto& pOrg : m_pOrganelles)
    {
        if (auto pMito = std::dynamic_pointer_cast<Mitochondrion>(pOrg))
        {
            return pMito;
        }
    }
    return nullptr;
}

bool Cell::consumeATP(double amount)
{
    // Consume ATP from medium at cell's position (center)
    float3 position(0.0f, 0.0f, 0.0f);
    return m_pMedium->consumeATP(amount, position);
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
            // Check both ATP and protein levels for transition
            if (cdk1 > 1000 && cyclinB > 1000 && consumeATP(ATPCosts::CHROMOSOME_CONDENSATION))
            {
                m_cellCycleState = CellCycleState::PROPHASE;
            }
            break;

        case CellCycleState::PROPHASE:
            // Transition to metaphase requires energy for spindle formation
            if (consumeATP(ATPCosts::SPINDLE_FORMATION))
            {
                // TODO: Add chromosome state monitoring
                // m_cellCycleState = CellCycleState::METAPHASE;
            }
            break;

        case CellCycleState::METAPHASE:
            // Transition to anaphase requires initial energy for chromosome movement
            if (consumeATP(ATPCosts::CHROMOSOME_MOVEMENT))
            {
                // TODO: Add spindle checkpoint monitoring
                // m_cellCycleState = CellCycleState::ANAPHASE;
            }
            break;

        case CellCycleState::ANAPHASE:
            // Continuous ATP consumption for chromosome movement
            if (consumeATP(ATPCosts::CHROMOSOME_MOVEMENT))
            {
                // TODO: Add chromosome position monitoring
                // m_cellCycleState = CellCycleState::TELOPHASE;
            }
            break;

        case CellCycleState::TELOPHASE:
            // Nuclear envelope reformation requires membrane fusion energy
            if (consumeATP(ATPCosts::MEMBRANE_FUSION))
            {
                // TODO: Add nuclear envelope monitoring
                // m_cellCycleState = CellCycleState::CYTOKINESIS;
            }
            break;

        case CellCycleState::CYTOKINESIS:
            // Cell membrane division requires fusion energy
            if (consumeATP(ATPCosts::MEMBRANE_FUSION))
            {
                // TODO: Add cytokinesis completion monitoring
                // m_cellCycleState = CellCycleState::INTERPHASE;
            }
            break;
    }
}
