#include "pch.h"
#include "Cell.h"
#include "Nucleus.h"
#include "Mitochondrion.h"
#include "Spindle.h"
#include "MRNA.h"
#include "log/ILog.h"

Cell::Cell(std::shared_ptr<Medium> pMedium, CellType type, size_t numChromosomes)
    : m_pMedium(pMedium)
    , m_cellCycleState(CellCycleState::INTERPHASE)
    , m_type(type)
{
    // Create DNA for the nucleus
    auto pDNA = std::make_shared<DNA>();
    
    // Create organelles
    m_pOrganelles.push_back(std::make_shared<Nucleus>(pDNA, numChromosomes));
    m_pOrganelles.push_back(std::make_shared<Mitochondrion>());
    // add other organelles as needed
}

void Cell::update(double fDt)
{
    // Update all organelles
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(fDt, *this, *m_pMedium);
    }
    
    // Check for cell cycle transitions based on conditions
    checkCellCycleTransitions();
    
    // Update medium
    m_pMedium->update(fDt);
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

bool Cell::consumeATP(double fAmount)
{
    // Consume ATP from medium at cell's position (center)
    float3 position(0.0f, 0.0f, 0.0f);
    return m_pMedium->consumeATP(fAmount, position);
}

void Cell::checkCellCycleTransitions()
{
    // Get key protein concentrations
    float3 center(0, 0, 0);
    double fCdk1 = m_pMedium->getProteinNumber("CDK-1", center);
    double fCyclinB = m_pMedium->getProteinNumber("CYB-1", center);
    double fPlk1 = m_pMedium->getProteinNumber("PLK-1", center);
    
    // Check conditions for each transition
    switch (m_cellCycleState)
    {
        case CellCycleState::INTERPHASE:
            // Check both ATP and protein levels for transition
            if (fCdk1 > 1000 && fCyclinB > 1000 && consumeATP(ATPCosts::fCHROMOSOME_CONDENSATION))
            {
                LOG_INFO("Cell switches from INTERPHASE to PROPHASE");
                m_cellCycleState = CellCycleState::PROPHASE;
                createSpindle();  // Create spindle as we enter prophase
            }
            break;

        case CellCycleState::PROPHASE:
            // Transition to metaphase requires energy for spindle formation
            if (consumeATP(ATPCosts::fSPINDLE_FORMATION))
            {
                if (auto pSpindle = getSpindle())
                {
                    if (pSpindle->isAssembled())
                    {
                        LOG_INFO("Cell switches from PROPHASE to METAPHASE");
                        m_cellCycleState = CellCycleState::METAPHASE;
                    }
                }
            }
            break;

        case CellCycleState::METAPHASE:
            // Transition to anaphase requires initial energy for chromosome movement
            if (consumeATP(ATPCosts::fCHROMOSOME_MOVEMENT))
            {
                // TODO: Add spindle checkpoint monitoring
                LOG_INFO("Cell switches from METAPHASE to ANAPHASE");
                m_cellCycleState = CellCycleState::ANAPHASE;
            }
            break;

        case CellCycleState::ANAPHASE:
            // Continuous ATP consumption for chromosome movement
            if (consumeATP(ATPCosts::fCHROMOSOME_MOVEMENT))
            {
                // TODO: Add chromosome position monitoring
                LOG_INFO("Cell switches from ANAPHASE to TELOPHASE");
                m_cellCycleState = CellCycleState::TELOPHASE;
            }
            break;

        case CellCycleState::TELOPHASE:
            // Nuclear envelope reformation requires membrane fusion energy
            if (consumeATP(ATPCosts::fMEMBRANE_FUSION))
            {
                LOG_INFO("Cell switches from TELOPHASE to CYTOKINESIS");
                m_cellCycleState = CellCycleState::CYTOKINESIS;
            }
            break;

        case CellCycleState::CYTOKINESIS:
            // Cell membrane division requires fusion energy
            if (consumeATP(ATPCosts::fMEMBRANE_FUSION))
            {
                destroySpindle();  // Destroy spindle as we complete division
                LOG_INFO("Cell switches from CYTOKINESIS to INTERPHASE");
                m_cellCycleState = CellCycleState::INTERPHASE;
            }
            break;
    }
}

void Cell::createSpindle()
{
    // Only create if we don't already have one
    if (!getSpindle())
    {
        m_pOrganelles.push_back(std::make_shared<Spindle>(m_type));
    }
}

void Cell::destroySpindle()
{
    m_pOrganelles.erase(
        std::remove_if(m_pOrganelles.begin(), m_pOrganelles.end(),
            [](const std::shared_ptr<Organelle>& pOrg) {
                return std::dynamic_pointer_cast<Spindle>(pOrg) != nullptr;
            }),
        m_pOrganelles.end());
}

std::shared_ptr<Spindle> Cell::getSpindle() const
{
    for (const auto& pOrg : m_pOrganelles)
    {
        if (auto pSpindle = std::dynamic_pointer_cast<Spindle>(pOrg))
        {
            return pSpindle;
        }
    }
    return nullptr;
}
