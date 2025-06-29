#include "pch.h"
#include "Cell.h"
#include "Nucleus.h"
#include "Mitochondrion.h"
#include "Spindle.h"
#include "EReticulum.h"
#include "Centrosome.h"
#include "Cortex.h"
#include "chemistry/MRNA.h"
#include "utils/log/ILog.h"

std::shared_ptr<Cell> Cell::createCell(std::shared_ptr<Medium> pInternalMedium,
                                     const std::vector<Chromosome>& chromosomes, 
                                     CellType type)
{
    auto cell = std::shared_ptr<Cell>(new Cell(pInternalMedium, chromosomes, type));
    cell->initializeOrganelles();
    cell->initializeCortex();
    return cell;
}

Cell::Cell(std::shared_ptr<Medium> pInternalMedium, const std::vector<Chromosome>& chromosomes, CellType type)
    : m_pInternalMedium(pInternalMedium)
    , m_cellCycleState(CellCycleState::INTERPHASE)
    , m_type(type)
    , m_chromosomes(chromosomes)
    , m_pOrganelles(static_cast<size_t>(StringDict::ID::ORGANELLE_END) - static_cast<size_t>(StringDict::ID::ORGANELLE_START))
{
    // Organelles will be initialized in initializeOrganelles() and initializeCortex()
}

void Cell::initializeOrganelles()
{
    // Create organelles using vector indexing
    m_pOrganelles[getOrganelleIndex(StringDict::ID::ORGANELLE_NUCLEUS)] = 
        std::make_shared<Nucleus>(std::weak_ptr<Cell>(shared_from_this()), m_chromosomes);
    m_pOrganelles[getOrganelleIndex(StringDict::ID::ORGANELLE_MITOCHONDRION)] = 
        std::make_shared<Mitochondrion>(std::weak_ptr<Cell>(shared_from_this()));
    m_pOrganelles[getOrganelleIndex(StringDict::ID::ORGANELLE_ENDOPLASMIC_RETICULUM)] = 
        std::make_shared<EReticulum>(std::weak_ptr<Cell>(shared_from_this()));
    // add other organelles as needed
}

void Cell::initializeCortex()
{
    // Create cortex as an organelle
    m_pOrganelles[getOrganelleIndex(StringDict::ID::ORGANELLE_CORTEX)] = 
        std::make_shared<Cortex>(std::weak_ptr<Cell>(shared_from_this()));
        
    // Initialize binding sites in the cortex
    auto pCortex = getCortex();
    if (pCortex)
    {
        pCortex->initializeBindingSites(4000000.0);
    }
}

std::shared_ptr<Cortex> Cell::getCortex() const
{
    size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_CORTEX);
    if (index < m_pOrganelles.size() && m_pOrganelles[index])
    {
        return std::dynamic_pointer_cast<Cortex>(m_pOrganelles[index]);
    }
    return nullptr;
}

void Cell::update(double fDt)
{
    // Update all organelles (including cortex)
    for (auto& pOrg : m_pOrganelles)
    {
        if (pOrg) {
            pOrg->update(fDt, *this);
        }
    }
    
    // Check for cell cycle transitions based on conditions
    checkCellCycleTransitions();
}

std::shared_ptr<Mitochondrion> Cell::getMitochondrion() const
{
    size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_MITOCHONDRION);
    if (index < m_pOrganelles.size() && m_pOrganelles[index])
    {
        return std::dynamic_pointer_cast<Mitochondrion>(m_pOrganelles[index]);
    }
    return nullptr;
}

bool Cell::consumeATP(double fAmount)
{
    // Consume ATP from internal medium at cell's position (center)
    float3 position(0.0f, 0.0f, 0.0f);
    return m_pInternalMedium->consumeATP(fAmount, position);
}

void Cell::checkCellCycleTransitions()
{
    // Get key protein concentrations from internal medium
    float3 center(0, 0, 0);
    Medium& internalMedium = *m_pInternalMedium;
    double fCdk1 = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CDK_1), center);
    double fCyclinB = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CYB_1), center);
    double fPlk1 = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::PLK_1), center);
    
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
        size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_SPINDLE);
        m_pOrganelles[index] = std::make_shared<Spindle>(std::weak_ptr<Cell>(shared_from_this()), m_type);
    }
}

void Cell::destroySpindle()
{
    size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_SPINDLE);
    if (index < m_pOrganelles.size())
    {
        m_pOrganelles[index].reset();
    }
}

std::shared_ptr<Spindle> Cell::getSpindle() const
{
    size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_SPINDLE);
    if (index < m_pOrganelles.size() && m_pOrganelles[index])
    {
        return std::dynamic_pointer_cast<Spindle>(m_pOrganelles[index]);
    }
    return nullptr;
}

void Cell::addOrganelle(StringDict::ID id, std::shared_ptr<Organelle> pOrganelle)
{
    if (pOrganelle) {
        size_t index = getOrganelleIndex(id);
        if (index < m_pOrganelles.size()) {
            m_pOrganelles[index] = pOrganelle;
        }
    }
}

std::shared_ptr<Centrosome> Cell::getCentrosome() const
{
    size_t index = getOrganelleIndex(StringDict::ID::ORGANELLE_CENTROSOME);
    if (index < m_pOrganelles.size() && m_pOrganelles[index])
    {
        return std::dynamic_pointer_cast<Centrosome>(m_pOrganelles[index]);
    }
    return nullptr;
}
