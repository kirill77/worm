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
    
    // Generate random normalized direction
    m_vDir = float3(dis(gen), dis(gen), dis(gen));
    float length = sqrtf(m_vDir.x * m_vDir.x + m_vDir.y * m_vDir.y + m_vDir.z * m_vDir.z);
    if (length > 0.0f) {
        m_vDir = float3(m_vDir.x / length, m_vDir.y / length, m_vDir.z / length);
    } else {
        m_vDir = float3(1.0f, 0.0f, 0.0f); // Default direction
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
    
    // Convert spherical to Cartesian coordinates
    m_vPosMicroM = float3(
        r * std::sin(phi) * std::cos(theta),
        r * std::sin(phi) * std::sin(theta), 
        r * std::cos(phi)
    );
}

void Y_TuRC::update(double dtSec, const float3& centrosomeWorldPos, const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium)
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
    float pCatFree = static_cast<float>(MoleculeConstants::MT_CATASTROPHE_BASE_FREE);   // s^-1 baseline catastrophe (free tip)
    float pRes = 0.05f;              // s^-1 baseline rescue

    if (m_mtLengthMicroM > 0.0f)
    {
        // Sample concentrations at plus-end position
        float3 originWorld = centrosomeWorldPos + m_vPosMicroM;
        float3 tipWorld = originWorld + m_vDir * m_mtLengthMicroM;
        assert(pCortex && "pCortex must not be null in Y_TuRC::update");
        float3 tipNorm = pCortex->worldToNormalized(tipWorld);
        double tubAlpha = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, m_species), tipNorm);
        double tubBeta  = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, m_species), tipNorm);
        float tubDimer = static_cast<float>(std::min(tubAlpha, tubBeta));
        double air1 = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, m_species), tipNorm);
        // Cortex-contact increases catastrophe probability
        if (m_mtContactCortex)
            pCatFree = static_cast<float>(MoleculeConstants::MT_CATASTROPHE_BASE_CONTACT); // stronger catastrophic bias at contact

        // Local soluble tubulin coupling (Michaelis–Menten-like)
        const float K_tub = static_cast<float>(MoleculeConstants::MT_K_TUB);
        float vGrow = vGrowMax * (tubDimer / (K_tub + tubDimer));
        // Catastrophe/rescue modulation by tubulin and AIR-1
        const float K_air = static_cast<float>(MoleculeConstants::MT_K_AIR);
        float f_air_cat = 1.0f / (1.0f + static_cast<float>(air1) / K_air);
        float f_tub_cat = (K_tub / (K_tub + std::max(tubDimer, 0.0f)));
        float f_tub_res = 1.0f + (std::max(tubDimer, 0.0f) / (K_tub + std::max(tubDimer, 0.0f)));
        pCatFree *= f_air_cat * f_tub_cat;
        pRes *= f_tub_res;

        if (m_mtState == MTState::Growing)
        {
            m_mtLengthMicroM += static_cast<float>(vGrow * dtSec);
            // GTP-cap grows with addition
            m_mtCapLengthMicroM += static_cast<float>(vGrow * dtSec);
            // Enhanced cortex intersection with cached triangle information
            Cortex::CortexRay intersection(originWorld, m_vDir);
            if (pCortex->findClosestIntersection(intersection) && m_mtLengthMicroM > intersection.distance)
            {
                m_mtLengthMicroM = intersection.distance;
                if (!m_mtContactCortex) {
                    // First contact - attempt cortical binding directly
                    m_mtContactCortex = true;
                    attemptCorticalBindingWithIntersection(intersection, pCortex, internalMedium, dtSec, rand01);
                }
            }
            else {
                m_mtContactCortex = false;
            }
            // Hydrolyze cap behind the tip
            m_mtCapLengthMicroM = std::max(0.0f, m_mtCapLengthMicroM - static_cast<float>(kHyd * dtSec));
            // Catastrophe probability rises when cap is depleted
            float pCat = pCatFree;
            if (m_mtCapLengthMicroM < static_cast<float>(MoleculeConstants::MT_CAP_DEPLETED_THRESHOLD_MICROM))
                pCat *= static_cast<float>(MoleculeConstants::MT_CAP_DEPLETION_CATASTROPHE_MULT);
            if (rand01() < pCat * dtSec)
            {
                m_mtState = MTState::Shrinking;
            }
        }
        else if (m_mtState == MTState::Shrinking)
        {
            m_mtLengthMicroM -= static_cast<float>(vShrink * dtSec);
            // Cap collapses in shrink
            m_mtCapLengthMicroM = std::max(0.0f, m_mtCapLengthMicroM - static_cast<float>(kHyd * dtSec));
            if (m_mtLengthMicroM <= 0.0f)
            {
                m_mtLengthMicroM = 0.0f;
                m_mtRefractorySec = 1.0f; // wait before re-nucleation
                m_mtContactCortex = false;
                
                // If was bound, reset binding state (cortex cleanup is lazy)
                if (m_bindingStrength > 0.0) {
                    m_bindingStrength = 0.0;
                    m_bindingTime = 0.0;
                }
            }
            else if (rand01() < pRes * dtSec)
            {
                m_mtState = MTState::Growing;
            }
        }
        else if (m_mtState == MTState::Bound)
        {
            // Update binding time and check for unbinding
            m_bindingTime += dtSec;
            if (shouldUnbind(dtSec, rand01)) {
                // Unbind from cortex (cortex cleanup is lazy) 
                unbindFromCortex();
            }
        }
    }
    else
    {
        if (m_mtRefractorySec > 0.0f)
        {
            m_mtRefractorySec = std::max(0.0f, m_mtRefractorySec - static_cast<float>(dtSec));
        }
        else
        {
            // Nucleate a new MT
            m_mtState = MTState::Growing;
            m_mtLengthMicroM = 0.02f;
            m_mtCapLengthMicroM = 0.02f;
            m_mtContactCortex = false;
            
            // Reset binding parameters
            m_bindingStrength = 0.0;
            m_bindingTime = 0.0;
        }
    }
}

void Y_TuRC::bindToCortex(double dyneinConc)
{
    m_mtState = MTState::Bound;
    m_bindingStrength = dyneinConc;
    m_bindingTime = 0.0;
}

void Y_TuRC::unbindFromCortex()
{
    m_mtState = MTState::Shrinking;  // Unbinding typically triggers catastrophe
    m_mtContactCortex = false;
    m_bindingStrength = 0.0;
    m_bindingTime = 0.0;
}

bool Y_TuRC::shouldUnbind(double dtSec, const std::function<float()>& rand01) const
{
    if (m_mtState != MTState::Bound) return false;
    
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

void Y_TuRC::attemptCorticalBindingWithIntersection(const Cortex::CortexRay& intersection,
                                                   const std::shared_ptr<Cortex>& pCortex, Medium& internalMedium, 
                                                   double dtSec, const std::function<float()>& rand01)
{
    
    // Sample NMY-2 (myosin II) concentration at contact point for binding probability
    const bool isOnCortex = true;
    float3 contactNorm = pCortex->worldToNormalized(intersection.worldHitPoint, isOnCortex);
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

