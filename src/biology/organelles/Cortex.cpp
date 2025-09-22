#include "pch.h"
#include "Cortex.h"
#include "Medium.h"
#include "Cell.h"
#include "utils/log/ILog.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "physics/tensionSphere/TensionSphere.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>

Cortex::Cortex(std::weak_ptr<Cell> pCell, double fThickness)
    : Organelle(pCell)
    , m_fThickness(fThickness)
{
    // Set the binding surface type for cortex
    m_surfaceType = StringDict::ID::ORGANELLE_CORTEX;

    // Initialize the tension sphere using the current cell volume
    auto pOwnedCell = getCell();
    if (pOwnedCell)
    {
        double volume = pOwnedCell->getInternalMedium().getVolumeMicroM();
        m_pTensionSphere = std::make_shared<TensionSphere>(2, volume);
    }
}

void Cortex::update(double fDtSec, Cell& cell)
{
    // Update internal medium - its dynamics are independent of external medium
    cell.getInternalMedium().update(fDtSec);
    
    // Maintain and advance tension sphere / BVH state for the cortex
    // Assume tension sphere exists; update volume and advance simulation
    double volume = cell.getInternalMedium().getVolumeMicroM();
    m_pTensionSphere->setVolume(volume);
    m_pTensionSphere->makeTimeStep(fDtSec);

    if (!m_pCortexBVH)
    {
        auto pCortexMesh = m_pTensionSphere->getEdgeMesh();
        if (pCortexMesh)
        {
            m_pCortexBVH = std::make_shared<BVHMesh>(pCortexMesh);
            cell.setCortexBVH(m_pCortexBVH);
        }
    }

    // Note: In a more advanced implementation, this method could include:
    // - Membrane fluidity changes
    // - Lipid raft formation/movement
    // - Membrane protein reorganization
    // - Passive transport based on concentration gradients
    // - Signal transduction
}

bool Cortex::initializeBindingSites(double totalAmount)
{
    // Ensure tension sphere and its mesh exist
    if (!m_pTensionSphere) {
        LOG_ERROR("Cannot initialize binding sites: tension sphere is not initialized");
        return false;
    }

    auto pMesh = m_pTensionSphere->getEdgeMesh();
    if (!pMesh) {
        LOG_ERROR("Cannot initialize binding sites: cortex mesh is not available");
        return false;
    }

    const uint32_t triangleCount = pMesh->getTriangleCount();
    if (triangleCount == 0) {
        LOG_ERROR("Cannot initialize binding sites: cortex mesh has no triangles");
        return false;
    }

    // Build area-weighted CDF over triangles for uniform surface sampling
    std::vector<double> triangleAreas(triangleCount);
    for (uint32_t t = 0; t < triangleCount; ++t) {
        triangleAreas[t] = pMesh->calculateTriangleArea(t);
    }

    const double totalArea = std::accumulate(triangleAreas.begin(), triangleAreas.end(), 0.0);
    if (totalArea <= 0.0) {
        LOG_ERROR("Cannot initialize binding sites: total cortex surface area is non-positive");
        return false;
    }

    std::vector<double> cdf(triangleCount);
    double cumulative = 0.0;
    for (uint32_t t = 0; t < triangleCount; ++t) {
        cumulative += triangleAreas[t] / totalArea;
        cdf[t] = cumulative;
    }
    // Guard against numerical issues
    cdf.back() = 1.0;

    // Prepare RNG
    static std::mt19937 rng(std::random_device{}());
    std::uniform_real_distribution<double> uni01(0.0, 1.0);

    // Target number of binding sites
    const uint32_t sitesPerAxis = 20;
    const uint32_t totalSites = sitesPerAxis * sitesPerAxis * sitesPerAxis; // 8000

    // Calculate amount per binding site to distribute totalAmount evenly
    const uint32_t totalPositions = totalSites;
    double amountPerPosition = (totalPositions > 0) ? (totalAmount / static_cast<double>(totalPositions)) : 0.0;

    // Molecule context
    auto pCell = getCell();
    Species species = pCell ? pCell->getSpecies() : Species::GENERIC;
    Molecule cortexBindingMol(StringDict::ID::ORGANELLE_CORTEX, ChemicalType::PROTEIN, species);

    m_pBindingSites.clear();
    m_pBindingSites.reserve(totalSites);

    for (uint32_t i = 0; i < totalSites; ++i) {
        // Sample triangle index by area-weighted CDF
        double u = uni01(rng);
        auto it = std::lower_bound(cdf.begin(), cdf.end(), u);
        uint32_t triIdx = static_cast<uint32_t>(std::distance(cdf.begin(), it));

        // Sample uniform barycentric coordinates over the triangle
        // Method: (1 - sqrt(r1), sqrt(r1)*(1 - r2), sqrt(r1)*r2)
        double r1 = uni01(rng);
        double r2 = uni01(rng);
        double sr1 = std::sqrt(r1);
        double b0 = 1.0 - sr1;
        double b1 = sr1 * (1.0 - r2);
        double b2 = sr1 * r2;

        BindingSite site;
        site.triangleIndex = triIdx;
        site.barycentric = float3(static_cast<float>(b0), static_cast<float>(b1), static_cast<float>(b2));
        // Initialize population on this binding site
        Population pop(amountPerPosition);
        pop.setBound(true);
        site.m_bsMolecules.emplace(cortexBindingMol, pop);
        m_pBindingSites.push_back(site);
    }

    return true;
}