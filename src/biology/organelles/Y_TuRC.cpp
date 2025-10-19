#include "Y_TuRC.h"
#include "Centrosome.h"
#include "Cortex.h"
#include "Cell.h"
#include "chemistry/molecules/StringDict.h"
#include <random>
#include "chemistry/molecules/Molecule.h"
#include <cmath>
#include <numbers>
#include <algorithm>
#include "chemistry/molecules/simConstants.h"

Y_TuRC::Y_TuRC(std::weak_ptr<Centrosome> pCentrosome)
    : m_pCentrosome(pCentrosome)
{
    // Initialize random direction for microtubule nucleation
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dis(-1.0f, 1.0f);
    
    // Generate random normalized direction for nucleation
    m_vTipDir = float3(dis(gen), dis(gen), dis(gen));
    float length = sqrtf(m_vTipDir.x * m_vTipDir.x + m_vTipDir.y * m_vTipDir.y + m_vTipDir.z * m_vTipDir.z);
    if (length > 0.0f) {
        m_vTipDir = float3(m_vTipDir.x / length, m_vTipDir.y / length, m_vTipDir.z / length);
    } else {
        m_vTipDir = float3(1.0f, 0.0f, 0.0f); // Default direction
    }
    
    // Initialize species from owning cell if available
    if (auto centrosome = pCentrosome.lock()) {
        if (auto cell = centrosome->getCell()) {
            m_species = cell->getSpecies();
        }
    }

    // Position randomly within PCM radius from the centrosome using spherical distribution
    float pcmRadius = 0.5f;  // Default radius
    if (auto centrosome = pCentrosome.lock()) {
        pcmRadius = centrosome->getPCMRadius();
    }
    
    // Generate uniform point within sphere
    std::uniform_real_distribution<float> unitDis(0.0f, 1.0f);
    std::uniform_real_distribution<float> thetaDis(0.0f, 2.0f * std::numbers::pi_v<float>);
    std::uniform_real_distribution<float> phiDis(-1.0f, 1.0f);
    
    float r = pcmRadius * std::cbrt(unitDis(gen));  // Cube root for uniform volume distribution
    float theta = thetaDis(gen);  // Azimuthal angle
    float phi = std::acos(phiDis(gen));  // Polar angle
    
    // Convert spherical to Cartesian coordinates and store as origin point
    float3 nucleationPos = float3(
        r * std::sin(phi) * std::cos(theta),
        r * std::sin(phi) * std::sin(theta), 
        r * std::cos(phi)
    );
    
    // Initialize with origin and a small initial segment
    getPoints().push_back(nucleationPos);
    getPoints().push_back(nucleationPos + m_vTipDir * 0.02f);
}

void Y_TuRC::update(double dtSec, const float3& centrosomeCellPos, const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium)
{
    // Defaults (embedded constants; not user-configurable)
    auto rand01 = []() {
        static thread_local std::mt19937 gen(std::random_device{}());
        static thread_local std::uniform_real_distribution<float> dis(0.0f, 1.0f);
        return dis(gen);
    };
    // Baseline rates at 20–22 °C
    const float vGrowMax = static_cast<float>(MoleculeConstants::MT_VGROW_MAX_UM_PER_S);    // µm/s (at high tubulin)
    const float vShrink  = static_cast<float>(MoleculeConstants::MT_VSHRINK_UM_PER_S);      // µm/s
    const float kHyd = static_cast<float>(MoleculeConstants::MT_CAP_HYDROLYSIS_RATE_S); // s^-1 cap hydrolysis
    float catastropheRate = static_cast<float>(MoleculeConstants::MT_CATASTROPHE_BASE_FREE);   // s^-1 baseline catastrophe (free tip)
    float rescueRate = 0.05f;              // s^-1 baseline rescue

    // Sample concentrations at plus-end position
    float3 tipCell = centrosomeCellPos + getTipPosition();
    assert(pCortex && "pCortex must not be null in Y_TuRC::update");
    float3 tipNorm = pCortex->cellToNormalized(tipCell);
    double tubAlpha = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, m_species), tipNorm);
    double tubBeta  = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, m_species), tipNorm);
    float tubDimer = static_cast<float>(std::min(tubAlpha, tubBeta));
    double air1 = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, m_species), tipNorm);
    // Cortex-contact increases catastrophe probability
    if (m_mtContactCortex)
        catastropheRate = static_cast<float>(MoleculeConstants::MT_CATASTROPHE_BASE_CONTACT); // stronger catastrophic bias at contact

    // Local soluble tubulin coupling (Michaelis–Menten-like)
    const float K_tub = static_cast<float>(MoleculeConstants::MT_K_TUB);
    float vGrow = vGrowMax * (tubDimer / (K_tub + tubDimer));
    // Catastrophe/rescue modulation by tubulin and AIR-1
    const float K_air = static_cast<float>(MoleculeConstants::MT_K_AIR);
    float f_air_cat = 1.0f / (1.0f + static_cast<float>(air1) / K_air);
    float f_tub_cat = (K_tub / (K_tub + std::max(tubDimer, 0.0f)));
    float f_tub_res = 1.0f + (std::max(tubDimer, 0.0f) / (K_tub + std::max(tubDimer, 0.0f)));
    catastropheRate *= f_air_cat * f_tub_cat;
    rescueRate *= f_tub_res;

    if (getState() == MTState::Growing)
    {
        float growthThisStep = static_cast<float>(vGrow * dtSec);
        const float segmentLength = static_cast<float>(MoleculeConstants::MT_SEGMENT_LENGTH_MICROM);
        
        // Grow the last segment in the tip direction
        float3 currentTip = getPoints().back();
        float3 newTip = currentTip + m_vTipDir * growthThisStep;
        getPoints().back() = newTip;
        
        // Add new full segments when last segment exceeds segment length
        float lastSegLen = getLastSegmentLength();
        while (lastSegLen >= segmentLength)
        {
            // Make the last segment exactly segmentLength
            float3 prevPoint = getPoints()[getPoints().size() - 2];
            getPoints().back() = prevPoint + m_vTipDir * segmentLength;
            
            // Add a new segment starting from current tip
            float3 segmentEnd = getPoints().back();
            float3 newSegmentEnd = segmentEnd + m_vTipDir * (lastSegLen - segmentLength);
            getPoints().push_back(newSegmentEnd);
            
            lastSegLen = getLastSegmentLength();
        }
        
        // Check cortex intersection for the tip segment only
        size_t lastIdx = getPoints().size() - 1;
        float3 segStart = centrosomeCellPos + getPoints()[lastIdx - 1];
        float3 segEnd = centrosomeCellPos + getPoints()[lastIdx];
        float3 segDir = segEnd - segStart;
        float segLen = sqrtf(segDir.x * segDir.x + segDir.y * segDir.y + segDir.z * segDir.z);
        
        m_mtContactCortex = false;
        if (segLen > 0.0f) {
            segDir = float3(segDir.x / segLen, segDir.y / segLen, segDir.z / segLen);
            
            Cortex::CortexRay intersection(segStart, segDir);
            if (pCortex->findClosestIntersection(intersection) && intersection.distance < segLen)
            {
                // Truncate MT at cortex contact
                float3 contactPoint = segStart + segDir * intersection.distance - centrosomeCellPos;
                getPoints()[lastIdx] = contactPoint;
                
                // First contact - attempt cortical binding directly
                m_mtContactCortex = true;
                attemptCorticalBindingWithIntersection(intersection, pCortex, internalMedium, dtSec, rand01);
            }
        }
        
        // Catastrophe probability rises when cap is depleted
        float capLength = getCapLengthMicroM();
        float effectiveCatastropheRate = catastropheRate;
        if (capLength < static_cast<float>(MoleculeConstants::MT_CAP_DEPLETED_THRESHOLD_MICROM))
            effectiveCatastropheRate *= static_cast<float>(MoleculeConstants::MT_CAP_DEPLETION_CATASTROPHE_MULT);
        if (rand01() < effectiveCatastropheRate * dtSec)
        {
            setState(MTState::Shrinking);
        }
    }
    else if (getState() == MTState::Shrinking)
    {
        float shrinkageThisStep = static_cast<float>(vShrink * dtSec);
        
        // Shrink the last segment
        float lastSegLen = getLastSegmentLength();
        if (shrinkageThisStep >= lastSegLen)
        {
            // Remove segments until shrinkage is consumed
            shrinkageThisStep -= lastSegLen;
            while (getPoints().size() > 2 && shrinkageThisStep > 0.0f)
            {
                getPoints().pop_back();
                float segLen = getLastSegmentLength();
                if (shrinkageThisStep >= segLen)
                {
                    shrinkageThisStep -= segLen;
                }
                else
                {
                    // Partially shrink this segment
                    float3 prevPoint = getPoints()[getPoints().size() - 2];
                    float3 currentTip = getPoints().back();
                    float3 segDir = currentTip - prevPoint;
                    float currentLen = sqrtf(segDir.x * segDir.x + segDir.y * segDir.y + segDir.z * segDir.z);
                    if (currentLen > 0.0f) {
                        segDir = float3(segDir.x / currentLen, segDir.y / currentLen, segDir.z / currentLen);
                        getPoints().back() = prevPoint + segDir * (currentLen - shrinkageThisStep);
                        m_vTipDir = segDir; // Update tip direction
                    }
                    shrinkageThisStep = 0.0f;
                }
            }
            
            // If down to minimum size (2 points), enter refractory period
            if (getPoints().size() == 2)
            {
                float remainingLen = getLastSegmentLength();
                if (remainingLen <= 0.02f || shrinkageThisStep > 0.0f)
                {
                    // Reset to small stub
                    float3 origin = getPoints()[0];
                    getPoints()[1] = origin + m_vTipDir * 0.02f;
                    m_mtRefractorySec = 1.0f; // wait before re-nucleation
                    m_mtContactCortex = false;
                    
                    // If was bound, reset binding state
                    if (m_bindingStrength > 0.0) {
                        m_bindingStrength = 0.0;
                        m_bindingTime = 0.0;
                    }
                }
            }
        }
        else
        {
            // Partially shrink the last segment
            float3 prevPoint = getPoints()[getPoints().size() - 2];
            float3 currentTip = getPoints().back();
            float3 segDir = currentTip - prevPoint;
            if (lastSegLen > 0.0f) {
                segDir = float3(segDir.x / lastSegLen, segDir.y / lastSegLen, segDir.z / lastSegLen);
                getPoints().back() = prevPoint + segDir * (lastSegLen - shrinkageThisStep);
            }
        }
        
        if (rand01() < rescueRate * dtSec)
        {
            setState(MTState::Growing);
        }
    }
    else if (getState() == MTState::Bound)
    {
        // Update binding time and check for unbinding
        m_bindingTime += dtSec;
        if (shouldUnbind(dtSec, rand01)) {
            // Unbind from cortex (cortex cleanup is lazy) 
            unbindFromCortex();
        }
    }
    
    // Handle refractory period
    if (m_mtRefractorySec > 0.0f)
    {
        m_mtRefractorySec = std::max(0.0f, m_mtRefractorySec - static_cast<float>(dtSec));
        if (m_mtRefractorySec == 0.0f)
        {
            // Re-nucleate
            setState(MTState::Growing);
            m_mtContactCortex = false;
            
            // Reset to small segment in tip direction
            float3 origin = getPoints()[0];
            getPoints()[1] = origin + m_vTipDir * 0.02f;
            
            // Reset binding parameters
            m_bindingStrength = 0.0;
            m_bindingTime = 0.0;
        }
    }
}

void Y_TuRC::bindToCortex(double dyneinConc)
{
    setState(MTState::Bound);
    m_bindingStrength = dyneinConc;
    m_bindingTime = 0.0;
}

void Y_TuRC::unbindFromCortex()
{
    setState(MTState::Shrinking);  // Unbinding typically triggers catastrophe
    m_mtContactCortex = false;
    m_bindingStrength = 0.0;
    m_bindingTime = 0.0;
}

bool Y_TuRC::shouldUnbind(double dtSec, const std::function<float()>& rand01) const
{
    if (getState() != MTState::Bound) return false;
    
    // Compute unbinding probability based on binding strength
    const double baseUnbindRate = MoleculeConstants::MT_CORTEX_UNBIND_BASE;
    const double weakUnbindRate = MoleculeConstants::MT_CORTEX_UNBIND_WEAK;
    const double K_bind = MoleculeConstants::MT_CORTEX_K_BIND;
    
    // Stronger binding (higher dynein) = lower unbinding rate
    double unbindRate = baseUnbindRate;
    if (m_bindingStrength < K_bind) {
        // Weak binding - higher unbinding rate
        double weakFactor = 1.0 - (m_bindingStrength / K_bind);
        unbindRate = baseUnbindRate + weakUnbindRate * weakFactor;
    }
    
    return rand01() < unbindRate * dtSec;
}

float Y_TuRC::getCapLengthMicroM() const
{
    // Cap length is proportional to how much the MT has grown
    // For simplicity, cap length equals total MT length (grows with polymerization)
    // In reality, GTP cap hydrolyzes over time - this could be tracked separately
    // For now, we use a simple proxy: cap depletes when MT is short
    float mtLength = getMTLengthMicroM();
    const float kHyd = static_cast<float>(MoleculeConstants::MT_CAP_HYDROLYSIS_RATE_S);
    // Simple approximation: cap is proportional to recent growth
    // This is a placeholder - proper implementation would track cap separately
    return mtLength * 0.1f; // Rough proxy: 10% of MT length is capped
}

void Y_TuRC::attemptCorticalBindingWithIntersection(const Cortex::CortexRay& intersection,
                                                   const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium, 
                                                   double dtSec, const std::function<float()>& rand01)
{
    
    // Sample NMY-2 (myosin II) concentration at contact point for binding probability
    const bool isOnCortex = true;
    float3 contactNorm = pCortex->cellToNormalized(intersection.worldHitPoint, isOnCortex);
    double dyneinConc = internalMedium.getMoleculeConcentration(
        Molecule(StringDict::ID::NMY_2, ChemicalType::PROTEIN, m_species), contactNorm);
    
    // Michaelis-Menten-like binding probability based on motor concentration
    const double K_bind = MoleculeConstants::MT_CORTEX_K_BIND;
    const double bindRate = MoleculeConstants::MT_CORTEX_BIND_RATE;
    double bindingProb = (dyneinConc / (K_bind + dyneinConc)) * bindRate * dtSec;
    
    if (rand01() < bindingProb) {
        // Successfully bind to cortex using cached intersection data
        bindToCortex(dyneinConc);
    }
}

