#include "pch.h"
#include "Centrosome.h"
#include "Cell.h"
#include "Medium.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "utils/log/ILog.h"
#include "Y_TuRC.h"
#include "Cortex.h"

Centrosome::Centrosome(std::weak_ptr<Cell> pCell, const float3& vNormalizedPos)
    : Organelle(pCell)
    , m_mToParent(affine3::identity())  // Initialize as identity transform
    , m_isDuplicated(false)
    , m_duplicationTime(0.0)
    , m_fPCMRadiusMicroM(0.5f)  // Default PCM radius of 0.5 micrometers
{
    // Set the position component of the transform
    m_mToParent.m_translation = vNormalizedPos;
    

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
        // TODO: Duplication gating will be implemented later.
        // We currently focus on the first ~100 seconds, during which duplication should not occur.
        // Leave this block intentionally inactive for now.
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
    // Local γ-tubulin concentration (molecules per µm^3)
    double gammaConc = internalMedium.getMoleculeConcentration(
        Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, species),
        m_mToParent.m_translation);
    
    // γ-tubulin recruitment driven by PCM maturation and local pool
    // Work purely with concentrations; additions to the medium remain count-based but should be derived from reaction models elsewhere.
    // Maintain a bound concentration proxy that accumulates with PCM maturation and decays slowly
    const double k_rec = 0.15;   // per second (tunable)
    const double k_loss = 0.005; // per second (tunable)
    const double dGammaBound = (k_rec * m_pcmMaturation * gammaConc - k_loss * m_gammaBoundConc) * dt;
    m_gammaBoundConc = std::max(0.0, m_gammaBoundConc + dGammaBound);

    // Target: proportional to local γ-tubulin concentration and PCM maturation
    // Keep target dimensionless using a concentration sensitivity factor (empirical tuning)
    // Provide a small basal count when PCM is present to ensure non-zero TuRCs pre-duplication
    const double beta = 50.0; // tunable sensitivity to bound concentration
    const int basal = (m_pcmMaturation > 0.05) ? 1 : 0; // minimal TuRCs when PCM emerges
    int targetRingComplexes = static_cast<int>(std::round(std::max(0.0, m_gammaBoundConc) * beta * m_pcmMaturation)) + basal;
    if (targetRingComplexes < 0) targetRingComplexes = 0;
    if (targetRingComplexes < 0) targetRingComplexes = 0;
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

    // Step microtubule dynamics for each existing Y_TuRC
    // Compute centrosome world position once
    float3 centrosomeWorldPos = float3(0,0,0);
    std::shared_ptr<Cortex> pCortex;
    if (auto pCellPtr = getCell())
    {
        pCortex = std::dynamic_pointer_cast<Cortex>(pCellPtr->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
        if (pCortex)
        {
            centrosomeWorldPos = pCortex->normalizedToWorld(getNormalizedPosition());
        }
    }
    for (auto& pRing : m_pRingComplexes)
    {
        if (pRing)
        {
            pRing->update(dt, centrosomeWorldPos, pCortex, internalMedium);
        }
    }
}

void Centrosome::updatePCMMaturation(double dt, const Cell& cell, Medium& internalMedium)
{
    Species species = cell.getSpecies();
    auto localCount = [&](StringDict::ID id) {
        return internalMedium.getMoleculeConcentration(Molecule(id, ChemicalType::PROTEIN, species), m_mToParent.m_translation);
    };
    const double spd2 = localCount(StringDict::ID::SPD_2);
    const double spd5 = localCount(StringDict::ID::SPD_5);
    const double plk1 = localCount(StringDict::ID::PLK_1);
    const double air1 = localCount(StringDict::ID::AIR_1);
    // Normalize by cell-average concentrations to make saturations unitless and grid-resolution-invariant
    // Use concentration API at the cell center (approximate whole-cell average for now)
    const float3 center(0.0f, 0.0f, 0.0f);
    const double spd2Avg = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::SPD_2, ChemicalType::PROTEIN, species), center);
    const double spd5Avg = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::SPD_5, ChemicalType::PROTEIN, species), center);
    const double plk1Avg = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::PLK_1, ChemicalType::PROTEIN, species), center);
    const double air1Avg = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, species), center);
    const double eps = 1e-12;
    const double nSPD2 = (spd2Avg > eps) ? (spd2 / spd2Avg) : 0.0;
    const double nSPD5 = (spd5Avg > eps) ? (spd5 / spd5Avg) : 0.0;
    const double nPLK1 = (plk1Avg > eps) ? (plk1 / plk1Avg) : 0.0;
    const double nAIR1 = (air1Avg > eps) ? (air1 / air1Avg) : 0.0;

    // Hill-like saturations on normalized signals (K=1 in normalized units)
    const double fSPD2 = nSPD2 / (1.0 + nSPD2);
    const double fSPD5 = nSPD5 / (1.0 + nSPD5);
    const double fKin = 1.0 + 0.5 * (nPLK1 / (1.0 + nPLK1)) + 0.3 * (nAIR1 / (1.0 + nAIR1));
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
