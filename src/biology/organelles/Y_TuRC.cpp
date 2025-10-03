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

Y_TuRC::Y_TuRC(std::weak_ptr<Centrosome> pCentrosome)
    : m_pCentrosome(pCentrosome)
    , m_nAlphaTubulins(0)
    , m_nBetaTubulins(0)
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
    const float vGrowMax = 0.45f;    // µm/s (at high tubulin)
    const float vShrink = 0.9f;      // µm/s
    const float kHyd = 0.3f;         // s^-1 cap hydrolysis
    float pCatFree = 0.12f;          // s^-1 baseline catastrophe (free tip)
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
            pCatFree = 0.6f; // stronger catastrophic bias at contact

        // Local soluble tubulin coupling (Michaelis–Menten-like)
        const float K_tub = 50.0f;
        float vGrow = vGrowMax * (tubDimer / (K_tub + tubDimer));
        // Catastrophe/rescue modulation by tubulin and AIR-1
        const float K_air = 10.0f;
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
            // Clamp by cortex distance along the MT direction in world space
            float maxLen = pCortex->distanceToCortex(originWorld, m_vDir);
            if (maxLen > 0.0f && m_mtLengthMicroM > maxLen)
            {
                m_mtLengthMicroM = maxLen;
                m_mtContactCortex = true;
            }
            else
                m_mtContactCortex = false;
            // Hydrolyze cap behind the tip
            m_mtCapLengthMicroM = std::max(0.0f, m_mtCapLengthMicroM - static_cast<float>(kHyd * dtSec));
            // Catastrophe probability rises when cap is depleted
            float pCat = pCatFree;
            if (m_mtCapLengthMicroM < 0.02f) pCat *= 3.0f;
            if (rand01() < pCat * dtSec)
            {
                m_mtState = MTState::Shrinking;
            }
        }
        else // Shrinking
        {
            m_mtLengthMicroM -= static_cast<float>(vShrink * dtSec);
            // Cap collapses in shrink
            m_mtCapLengthMicroM = std::max(0.0f, m_mtCapLengthMicroM - static_cast<float>(kHyd * dtSec));
            if (m_mtLengthMicroM <= 0.0f)
            {
                m_mtLengthMicroM = 0.0f;
                m_mtRefractorySec = 1.0f; // wait before re-nucleation
                m_mtContactCortex = false;
            }
            else if (rand01() < pRes * dtSec)
            {
                m_mtState = MTState::Growing;
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
        }
    }
}