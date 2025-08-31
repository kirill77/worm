#include "pch.h"
#include "Centrosome.h"
#include "Cell.h"
#include "Medium.h"
#include "chemistry/MoleculeWiki.h"
#include "utils/log/ILog.h"
#include "Y_TuRC.h"

Centrosome::Centrosome(std::weak_ptr<Cell> pCell, const float3& vNormalizedPos)
    : Organelle(pCell)
    , m_mToParent(affine3::identity())  // Initialize as identity transform
    , m_isDuplicated(false)
    , m_duplicationTime(0.0)
    , m_fPCMRadiusMicroM(0.5f)  // Default PCM radius of 0.5 micrometers
{
    // Set the position component of the transform
    m_mToParent.m_translation = vNormalizedPos;
    
    // Set the binding surface type to CENTROSOME
    m_surfaceType = StringDict::ID::ORGANELLE_CENTROSOME;
    
    // Initialize centrosome-specific proteins
    if (auto pCellPtr = pCell.lock()) {
        auto& internalMedium = pCellPtr->getInternalMedium();
        
        // γ-tubulin is now produced through gene expression, not constant addition
        
        // Add other centrosome proteins
        MPopulation pericentrin(StringDict::idToString(StringDict::ID::PERICENTRIN), 500.0);
        internalMedium.addMolecule(pericentrin, m_mToParent.m_translation);
        
        MPopulation ninein(StringDict::idToString(StringDict::ID::NINEIN), 300.0);
        internalMedium.addMolecule(ninein, m_mToParent.m_translation);
    }
    
    // Note: Ring complexes will be created in the update method
}

void Centrosome::update(double dt, Cell& cell)
{
    // Get the cell's internal medium
    auto& internalMedium = cell.getInternalMedium();
    
    // Check for centrosome duplication during S phase
    if (!m_isDuplicated) {
        // Check if we're in S phase (high CDK-2 and Cyclin E levels)
        double cdk2 = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CDK_2), m_mToParent.m_translation);
        double cyclinE = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CCE_1), m_mToParent.m_translation);
        
        // Trigger duplication when CDK-2 and Cyclin E are high enough
        if (cdk2 > 800.0 && cyclinE > 800.0) {
            duplicate();
            LOG_INFO("Centrosome duplication triggered at position (%.2f, %.2f, %.2f)", 
                     m_mToParent.m_translation.x, m_mToParent.m_translation.y, m_mToParent.m_translation.z);
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
                float3 polePosition = m_mToParent.m_translation;
                if (m_mToParent.m_translation.y > 0) {
                    polePosition.y = 0.8f;  // Anterior pole
                } else {
                    polePosition.y = -0.8f; // Posterior pole
                }
                m_mToParent.m_translation = polePosition;
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
            break;
            
        default:
            // In interphase, centrosome stays near the nucleus
            break;
    }
    
    // γ-tubulin is now regulated through transcription, not constant addition
    
    // Manage ring complexes based on gamma-tubulin levels
    double gammaTubulinCount = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::GAMMA_TUBULIN), m_mToParent.m_translation);
    
    // Target: ~1 ring complex per 50 gamma-tubulin proteins
    int targetRingComplexes = static_cast<int>(gammaTubulinCount / 50.0);
    int currentRingComplexes = static_cast<int>(m_pRingComplexes.size());
    
    if (currentRingComplexes < targetRingComplexes) {
        // Create new ring complexes using shared_from_this() and cast to weak_ptr
        std::weak_ptr<Centrosome> thisWeakPtr = std::static_pointer_cast<Centrosome>(shared_from_this());
        for (int i = currentRingComplexes; i < targetRingComplexes; i++) {
            auto ringComplex = std::make_shared<Y_TuRC>(thisWeakPtr);
            m_pRingComplexes.push_back(ringComplex);
        }
    } else if (currentRingComplexes > targetRingComplexes) {
        // Remove excess ring complexes
        int complexesToRemove = currentRingComplexes - targetRingComplexes;
        for (int i = 0; i < complexesToRemove; i++) {
            m_pRingComplexes.pop_back();
        }
    }
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
            internalMedium.addMolecule(gammaTubulin, m_mToParent.m_translation);
            
            // Add other duplication-related proteins
            MPopulation plk4(StringDict::idToString(StringDict::ID::PLK_4), 200.0);  // Polo-like kinase 4 for centriole duplication
            internalMedium.addMolecule(plk4, m_mToParent.m_translation);
            
            // Add additional ring complexes for the duplicated centrosome
            int additionalRingComplexes = static_cast<int>(m_pRingComplexes.size() * 0.5); // 50% more
            std::weak_ptr<Centrosome> thisWeakPtr = std::static_pointer_cast<Centrosome>(shared_from_this());
            for (int i = 0; i < additionalRingComplexes; i++) {
                auto ringComplex = std::make_shared<Y_TuRC>(thisWeakPtr);
                m_pRingComplexes.push_back(ringComplex);
            }
            
            LOG_INFO("Added %d additional ring complexes during centrosome duplication", additionalRingComplexes);
        }
    }
}
