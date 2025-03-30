#include "pch.h"
#include "Cell.h"
#include "Nucleus.h"
#include "Mitochondrion.h"
#include "Spindle.h"
#include "MRNA.h"
#include "log/ILog.h"

Cell::Cell(std::shared_ptr<Cortex> pCortex, const std::vector<Chromosome>& chromosomes, CellType type)
    : m_pCortex(pCortex)
    , m_cellCycleState(CellCycleState::INTERPHASE)
    , m_type(type)
{
    // Create organelles
    m_pOrganelles.push_back(std::make_shared<Nucleus>(chromosomes));
    m_pOrganelles.push_back(std::make_shared<Mitochondrion>());
    // add other organelles as needed
    
    // Initialize binding sites in the cell's membrane
    if (m_pCortex)
    {
        m_pCortex->initializeBindingSites(4000000.0);
    }
}

void Cell::update(double fDt)
{
    // Update all organelles - pass the internal medium to organelles
    Medium& internalMedium = m_pCortex->getInternalMedium();
    for (auto& pOrg : m_pOrganelles)
    {
        pOrg->update(fDt, *this, internalMedium);
    }
    
    // Check for cell cycle transitions based on conditions
    checkCellCycleTransitions();
    
    // Update the membrane which in turn will update the internal medium
    m_pCortex->update(fDt);
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
    // Consume ATP from internal medium at cell's position (center)
    float3 position(0.0f, 0.0f, 0.0f);
    return m_pCortex->getInternalMedium().consumeATP(fAmount, position);
}

void Cell::checkCellCycleTransitions()
{
    // Get key protein concentrations from internal medium
    float3 center(0, 0, 0);
    Medium& internalMedium = m_pCortex->getInternalMedium();
    double fCdk1 = internalMedium.getProteinNumber("CDK-1", center);
    double fCyclinB = internalMedium.getProteinNumber("CYB-1", center);
    double fPlk1 = internalMedium.getProteinNumber("PLK-1", center);
    
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
