#include "pch.h"
#include "Centrosome.h"
#include "Cell.h"
#include "Medium.h"
#include "chemistry/molecules/MoleculeWiki.h"
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
        MPopulation pericentrin(Molecule(StringDict::ID::PERICENTRIN, ChemicalType::PROTEIN), 500.0);
        internalMedium.addMolecule(pericentrin, m_mToParent.m_translation);
        
        MPopulation ninein(Molecule(StringDict::ID::NINEIN, ChemicalType::PROTEIN), 300.0);
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
        double cdk2 = internalMedium.getMoleculeNumber(Molecule(StringDict::ID::CDK_2, ChemicalType::PROTEIN), m_mToParent.m_translation);
        double cyclinE = internalMedium.getMoleculeNumber(Molecule(StringDict::ID::CCE_1, ChemicalType::PROTEIN), m_mToParent.m_translation);
        
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

    updatePCMMaturation(dt, cell, internalMedium);

    updateGammaAndRingComplexes(dt, cell, internalMedium);
}

void Centrosome::updateGammaAndRingComplexes(double dt, const Cell& cell, Medium& internalMedium)
{
    Species species = cell.getSpecies();
    double gammaTubulinCount = internalMedium.getMoleculeNumber(
        Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, species),
        m_mToParent.m_translation);

    // γ-tubulin recruitment driven by PCM maturation and local cytosolic pool
    const double gammaCyt = gammaTubulinCount; // proxy for local available pool; later separate bound/unbound
    const double k_rec = 0.1;   // per second (tunable)
    const double k_loss = 0.01; // per second (tunable)
    const double dGamma = (k_rec * m_pcmMaturation * (gammaCyt) - k_loss * gammaTubulinCount) * dt;
    if (dGamma > 0.0)
    {
        MPopulation gammaAdd(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, species), dGamma);
        internalMedium.addMolecule(gammaAdd, m_mToParent.m_translation);
        gammaTubulinCount += dGamma;
    }

    // Target: proportional to γ-tubulin and PCM maturation
    const double tuRCPerGamma = 1.0 / 50.0;
    int targetRingComplexes = static_cast<int>(gammaTubulinCount * tuRCPerGamma * m_pcmMaturation);
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

void Centrosome::updatePCMMaturation(double dt, const Cell& cell, Medium& internalMedium)
{
    Species species = cell.getSpecies();
    auto localCount = [&](StringDict::ID id) {
        return internalMedium.getMoleculeNumber(Molecule(id, ChemicalType::PROTEIN, species), m_mToParent.m_translation);
    };
    const double spd2 = localCount(StringDict::ID::SPD_2);
    const double spd5 = localCount(StringDict::ID::SPD_5);
    const double plk1 = localCount(StringDict::ID::PLK_1);
    const double air1 = localCount(StringDict::ID::AIR_1);
    // Hill-like saturations (tunable)
    const double K2 = 200.0, K5 = 200.0, Kp = 100.0, Ka = 50.0;
    const double fSPD2 = spd2 / (K2 + std::max(1.0, spd2));
    const double fSPD5 = spd5 / (K5 + std::max(1.0, spd5));
    const double fKin = 1.0 + 0.5 * (plk1 / (Kp + std::max(1.0, plk1))) + 0.3 * (air1 / (Ka + std::max(1.0, air1)));
    // Kinetics (tunable)
    const double k_on = 0.2;   // per second
    const double k_off = 0.02; // per second
    const double pcmProd = k_on * fSPD2 * fSPD5 * fKin * (1.0 - m_pcmMaturation);
    const double pcmLoss = k_off * m_pcmMaturation;
    m_pcmMaturation = std::clamp(m_pcmMaturation + (pcmProd - pcmLoss) * dt, 0.0, 1.0);
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
            MPopulation gammaTubulin(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN), 500.0);
            internalMedium.addMolecule(gammaTubulin, m_mToParent.m_translation);
            
            // Add other duplication-related proteins
            MPopulation plk4(Molecule(StringDict::ID::PLK_4, ChemicalType::PROTEIN), 200.0);  // Polo-like kinase 4 for centriole duplication
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
