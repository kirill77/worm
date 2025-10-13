#include "pch.h"
#include "Cortex.h"
#include "Medium.h"
#include "Cell.h"
#include "utils/log/ILog.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "geometry/geomHelpers/BVHCache.h"
#include "TensionSphere.h"
#include "geometry/BVH/BVH.h"
#include "geometry/BVH/ITraceableObject.h"
#include <algorithm>
#include <cmath>
#include <random>
#include <numeric>
#include <limits>

Cortex::Cortex(std::weak_ptr<Cell> pCell, double fThickness)
    : Organelle(pCell)
    , m_fThickness(fThickness)
{

    // Initialize the tension sphere using the current cell volume
    auto pOwnedCell = getCell();
    if (pOwnedCell)
    {
        double volume = pOwnedCell->getInternalMedium().getVolumeMicroM();
        m_pTensionSphere = std::make_shared<TensionSphere>(2, volume);
        // Initialize BVH from the tension sphere mesh for spatial queries
        if (m_pTensionSphere)
        {
            auto pMesh = m_pTensionSphere->getEdgeMesh();
            if (pMesh)
            {
                m_pCortexBVH = BVHCache::instance().getOrCreate(pMesh);

                // Validate mapping consistency between normalizedToWorld and worldToNormalized
                static std::mt19937 rng(std::random_device{}());
                std::uniform_real_distribution<float> uni(-1.0f, 1.0f);
                for (int i = 0; i < 10; ++i)
                {
                    float3 n = float3(uni(rng), uni(rng), uni(rng));
                    float r = length(n);
                    if (r < 1e-6f) n = float3(1, 0, 0), r = 1.0f;
                    float3 w = normalizedToWorld(n);
                    float3 n2 = worldToNormalized(w);
                    float3 d = n2 - n;
                    float err = length(d);
                    assert(err < 1e-3f && "Cortex world/normalized mappings must be approximately inverse after clamping to unit length");
                }
            }
        }
    }

    // Initialize list of molecules that can bind to cortex
    // For now, include the cortex-binding protein key
    Species species = pOwnedCell ? pOwnedCell->getSpecies() : Species::GENERIC;
    m_bindableMolecules.emplace_back(StringDict::ID::ORGANELLE_CORTEX, ChemicalType::PROTEIN, species);
    // Also include specific cortex-bound PAR complexes
    m_bindableMolecules.emplace_back(StringDict::ID::PAR_1_CORTEX, ChemicalType::PROTEIN, species);
    m_bindableMolecules.emplace_back(StringDict::ID::PAR_2_CORTEX, ChemicalType::PROTEIN, species);
    m_bindableMolecules.emplace_back(StringDict::ID::PAR_3_CORTEX, ChemicalType::PROTEIN, species);
}

void Cortex::update(double fDtSec, Cell& cell)
{
    // Pull molecules from grid to binding sites prior to shape update
    pullBindingSiteMoleculesFromMedium();

    // Maintain and advance tension sphere / BVH state for the cortex
    // Assume tension sphere exists; update volume and advance simulation
    double volume = cell.getInternalMedium().getVolumeMicroM();
    m_pTensionSphere->setVolume(volume);
    m_pTensionSphere->makeTimeStep(fDtSec);

    // Ensure BVHMesh is up-to-date with the latest mesh version
    auto pMesh = m_pTensionSphere->getEdgeMesh();
    m_pCortexBVH = BVHCache::instance().getOrCreate(pMesh);

    // Push molecules back into the grid at updated positions after shape update
    transferBindingSiteMoleculesToMedium();

    // After shape update, update per-cell volumes for concentration queries
    cell.getInternalMedium().updateGridCellVolumes(*this);

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

        CortexMolecules site;
        site.m_triangleIndex = triIdx;
        site.setBarycentric(float3(static_cast<float>(b0), static_cast<float>(b1), static_cast<float>(b2)));
        // Initialize population on this binding site
        Population pop(amountPerPosition);
        pop.setBound(true);
        site.m_bsMolecules.emplace(cortexBindingMol, pop);
        m_pBindingSites.push_back(site);
    }

    // After creating binding sites, move their molecules into the simulation grid
    transferBindingSiteMoleculesToMedium();

    return true;
}

void Cortex::transferBindingSiteMoleculesToMedium()
{
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot transfer binding site molecules: cell reference is invalid");
        return;
    }
    if (!m_pTensionSphere) {
        LOG_ERROR("Cannot transfer binding site molecules: tension sphere is not initialized");
        return;
    }

    auto pMesh = m_pTensionSphere->getEdgeMesh();
    if (!pMesh) {
        LOG_ERROR("Cannot transfer binding site molecules: cortex mesh is not available");
        return;
    }

    Medium& medium = pCell->getInternalMedium();

    const uint32_t triangleCount = pMesh->getTriangleCount();
    for (auto& site : m_pBindingSites) {
        if (site.m_triangleIndex >= triangleCount)
            continue;

        // Update normalized [-1,1] position on cortex from triangle + barycentric
        site.setNormalized(baryToNormalized(site.m_triangleIndex, site.getBarycentric()));
        const float3 pos = site.getNormalized();

        // Transfer each molecule population to the medium at this position and zero source immediately
        for (auto& kv : site.m_bsMolecules) {
            const Molecule& mol = kv.first;
            Population& pop = kv.second;
            if (pop.m_fNumber > 0.0) {
                MPopulation mpop(mol, pop);
                assert(mpop.isBound()); // we're working with molecules bound to cortex here
                medium.addMolecule(mpop, pos);
                // zero out source
                pop.m_fNumber = 0.0;
                pop.setBound(false);
            }
        }
    }
}

void Cortex::pullBindingSiteMoleculesFromMedium()
{
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot pull binding site molecules: cell reference is invalid");
        return;
    }

    // m_normalized should not change while molecules are in the medium.
    // Simply delegate moving molecules from grid cells into binding sites.
    Medium& medium = pCell->getInternalMedium();
    medium.toBindingSites(m_pBindingSites, m_bindableMolecules);
}



float3 Cortex::normalizedToWorld(const float3& normalizedPos)
{
    assert(m_pCortexBVH && "Cortex BVH must be initialized before normalizedToWorld");

    // Get bounding box center (for ray origin)
    const box3 bbox = m_pCortexBVH->getBox();
    const float3 center = bbox.center();

    // Early-out for near-origin input
    float nLen = length(normalizedPos);
    if (nLen < 1e-6f)
        return center;

    // Decompose normalized input into (unit direction, scalar s in [0,1])
    // s encodes the fraction to the cortex along the ray; use L-infinity norm to remain inside cube
    float s = std::max(std::max(std::abs(normalizedPos.x), std::abs(normalizedPos.y)), std::abs(normalizedPos.z));
    s = std::min(1.0f, std::max(0.0f, s));
    float3 dirInf = normalizedPos / nLen; // unit direction in normalized space (L2)

    // Convert normalized-space direction to world-space ray direction by applying half-extents, then renormalizing
    const float3 half = (bbox.m_maxs - bbox.m_mins) * 0.5f;
    float3 dirWorldPre = float3(dirInf.x * half.x, dirInf.y * half.y, dirInf.z * half.z);
    float preLen = length(dirWorldPre);
    if (preLen < 1e-6f)
        return center;
    float3 dirWorldUnit = dirWorldPre / preLen;

    // Distance to cortex along this direction
    CortexRay ray(center, dirWorldUnit);
    float distCortex = findClosestIntersection(ray) ? ray.getDistance() : 0.0f;
    if (distCortex <= 0.0f)
        return center; // Degenerate case

    // Map to world along the same ray by fraction s of the cortex distance
    return center + dirWorldUnit * (distCortex * s);
}


float3 Cortex::worldToNormalized(const float3& worldPos, bool isOnCortex) const
{
    assert(m_pCortexBVH && "Cortex BVH must be initialized before worldToNormalized");
    const box3 bbox = m_pCortexBVH->getBox();
    const float3 center = bbox.center();
    float3 v = worldPos - center;
    float len = length(v);
    if (len < 1e-6f) return float3(0, 0, 0);
    float3 dirWorldUnit = v / len;
    float distCortex = 0.0f;
    if (isOnCortex) {
        distCortex = len;
    }
    else {
        CortexRay ray(center, dirWorldUnit);
        distCortex = findClosestIntersection(ray) ? ray.getDistance() : 0.0f;
    }
    if (distCortex <= 0.0f)
        return float3(0, 0, 0);

    // Fraction along the ray to the cortex (clamped)
    float s = std::min(1.0f, std::max(0.0f, len / distCortex));

    // Undo anisotropic scaling to recover normalized-space direction, then L-infinity normalize
    const float3 half = (bbox.m_maxs - bbox.m_mins) * 0.5f;
    const float eps = 1e-8f;
    float3 pre = float3(
        (std::abs(half.x) > eps) ? (dirWorldUnit.x / half.x) : 0.0f,
        (std::abs(half.y) > eps) ? (dirWorldUnit.y / half.y) : 0.0f,
        (std::abs(half.z) > eps) ? (dirWorldUnit.z / half.z) : 0.0f
    );
    float maxAbs = std::max(std::abs(pre.x), std::max(std::abs(pre.y), std::abs(pre.z)));
    if (maxAbs < eps) return float3(0, 0, 0);
    float3 dirInf = pre / maxAbs;

    return dirInf * s;
}

float3 Cortex::baryToNormalized(uint32_t triangleIndex, const float3& barycentric) const
{
    if (!m_pTensionSphere)
        return float3(0, 0, 0);
    auto pMesh = m_pTensionSphere->getEdgeMesh();
    if (!pMesh)
        return float3(0, 0, 0);
    if (triangleIndex >= pMesh->getTriangleCount())
        return float3(0, 0, 0);

    const uint3 tri = pMesh->getTriangleVertices(triangleIndex);
    const float3 v0 = pMesh->getVertexPosition(tri.x);
    const float3 v1 = pMesh->getVertexPosition(tri.y);
    const float3 v2 = pMesh->getVertexPosition(tri.z);

    const float3 world = v0 * barycentric.x + v1 * barycentric.y + v2 * barycentric.z;
    const bool isOnCortex = true;
    return worldToNormalized(world, isOnCortex);
}

bool Cortex::findClosestIntersection(CortexRay& ray) const
{
    assert(m_pCortexBVH && "Cortex BVH must be initialized before findClosestIntersection");
    // Trace directly using the ray object (implements IRay)
    // BVHMesh is refreshed via BVHCache in update(); its BVH is up-to-date here
    const BVH& bvh = m_pCortexBVH->getBVH();
    bvh.trace(ray, 0);
    return ray.hasHit;
}

// Cortex::CortexRay definitions
Cortex::CortexRay::CortexRay(const float3& origin, const float3& direction)
{
    m_vPos = origin;
    m_vDir = normalize(direction);
    m_fMin = 0.0f;
    m_fMax = std::numeric_limits<float>::max();
}

void Cortex::CortexRay::notifyIntersection(float fDist, const ITraceableObject*, uint32_t uSubObj)
{
    if (fDist >= m_fMin && fDist <= m_fMax && fDist < distance) {
        distance = fDist;
        triangleIndex = uSubObj;
        hasHit = true;
        worldHitPoint = m_vPos + m_vDir * fDist;
    }
}

