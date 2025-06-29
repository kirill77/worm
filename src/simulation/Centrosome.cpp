#include "pch.h"
#include "Centrosome.h"
#include "Cell.h"
#include "Medium.h"
#include "molecules/ProteinWiki.h"
#include "utils/log/ILog.h"

Centrosome::Centrosome(std::weak_ptr<Cell> pCell, const float3& position)
    : Organelle(pCell)
    , m_position(position)
    , m_isDuplicated(false)
    , m_duplicationTime(0.0)
{
    // Set the binding surface type to CENTROSOME
    m_surfaceType = StringDict::ID::BS_CENTROSOME;
    
    // Initialize centrosome-specific proteins
    if (auto pCellPtr = pCell.lock()) {
        auto& internalMedium = pCellPtr->getInternalMedium();
        
        // Add centrosome-specific proteins like γ-tubulin
        MPopulation gammaTubulin(StringDict::idToString(StringDict::ID::GAMMA_TUBULIN), 1000.0);
        internalMedium.addProtein(gammaTubulin, m_position);
        
        // Add other centrosome proteins
        MPopulation pericentrin(StringDict::idToString(StringDict::ID::PERICENTRIN), 500.0);
        internalMedium.addProtein(pericentrin, m_position);
        
        MPopulation ninein(StringDict::idToString(StringDict::ID::NINEIN), 300.0);
        internalMedium.addProtein(ninein, m_position);
    }
}

void Centrosome::update(double dt, Cell& cell)
{
    // Get the cell's internal medium
    auto& internalMedium = cell.getInternalMedium();
    
    // Check for centrosome duplication during S phase
    if (!m_isDuplicated) {
        // Check if we're in S phase (high CDK-2 and Cyclin E levels)
        double cdk2 = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CDK_2), m_position);
        double cyclinE = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CCE_1), m_position);
        
        // Trigger duplication when CDK-2 and Cyclin E are high enough
        if (cdk2 > 800.0 && cyclinE > 800.0) {
            duplicate();
            LOG_INFO("Centrosome duplication triggered at position (%.2f, %.2f, %.2f)", 
                     m_position.x, m_position.y, m_position.z);
        }
    }
    
    // Update centrosome position based on cell cycle state
    CellCycleState cellCycleState = cell.getCellCycleState();
    
    switch (cellCycleState) {
        case CellCycleState::PROPHASE:
        case CellCycleState::METAPHASE:
            // During mitosis, centrosomes move to opposite poles
            if (m_isDuplicated) {
                // Move to opposite poles (simplified - in reality this is more complex)
                float3 polePosition = m_position;
                if (m_position.y > 0) {
                    polePosition.y = 0.8f;  // Anterior pole
                } else {
                    polePosition.y = -0.8f; // Posterior pole
                }
                m_position = polePosition;
            }
            break;
            
        case CellCycleState::ANAPHASE:
        case CellCycleState::TELOPHASE:
            // Continue moving to poles during anaphase/telophase
            break;
            
        case CellCycleState::CYTOKINESIS:
            // Reset duplication state after cell division
            m_isDuplicated = false;
            m_duplicationTime = 0.0;
            // Reset position to center for next cycle
            m_position = float3(0, 0, 0);
            break;
            
        default:
            // In interphase, centrosome stays near the nucleus
            break;
    }
    
    // Update protein concentrations at the centrosome position
    // This ensures centrosome proteins are properly localized
    MPopulation gammaTubulin(StringDict::idToString(StringDict::ID::GAMMA_TUBULIN), 1000.0);
    internalMedium.addProtein(gammaTubulin, m_position);
}

void Centrosome::duplicate()
{
    if (!m_isDuplicated) {
        m_isDuplicated = true;
        m_duplicationTime = 0.0;  // Reset timer
        
        // Add duplicated centrosome proteins
        if (auto pCellPtr = getCell()) {
            auto& internalMedium = pCellPtr->getInternalMedium();
            
            // Add additional γ-tubulin for the duplicated centrosome
            MPopulation gammaTubulin(StringDict::idToString(StringDict::ID::GAMMA_TUBULIN), 500.0);
            internalMedium.addProtein(gammaTubulin, m_position);
            
            // Add other duplication-related proteins
            MPopulation plk4(StringDict::idToString(StringDict::ID::PLK_4), 200.0);  // Polo-like kinase 4 for centriole duplication
            internalMedium.addProtein(plk4, m_position);
        }
    }
}
