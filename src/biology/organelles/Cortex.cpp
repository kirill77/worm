#include "pch.h"
#include "Cortex.h"
#include "Medium.h"
#include "Cell.h"
#include "utils/log/ILog.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "physics/tensionSphere/TensionSphere.h"
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
    // Set the binding surface type for cortex
    m_surfaceType = StringDict::ID::ORGANELLE_CORTEX;

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
                m_pCortexBVH = std::make_shared<BVHMesh>(pMesh);
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

    if (!m_pCortexBVH)
    {
        auto pCortexMesh = m_pTensionSphere->getEdgeMesh();
        if (pCortexMesh)
        {
            m_pCortexBVH = std::make_shared<BVHMesh>(pCortexMesh);
        }
    }

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

        BindingSite site;
        site.m_triangleIndex = triIdx;
        site.m_barycentric = float3(static_cast<float>(b0), static_cast<float>(b1), static_cast<float>(b2));
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
        site.m_normalized = baryToNormalized(site.m_triangleIndex, site.m_barycentric);
        const float3 pos = site.m_normalized;

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

namespace {
class CortexRay : public IRay
{
public:
    float m_closestDistance;
    bool m_hasIntersection;

    CortexRay(const float3& pos, const float3& dir)
        : m_closestDistance(std::numeric_limits<float>::max())
        , m_hasIntersection(false)
    {
        m_vPos = pos;
        m_vDir = normalize(dir);
        m_fMin = 0.0f;
        m_fMax = std::numeric_limits<float>::max();
    }

    void notifyIntersection(float fDist, const ITraceableObject* /*pObject*/, uint32_t /*uSubObj*/) override
    {
        if (fDist >= m_fMin && fDist <= m_fMax && fDist < m_closestDistance)
        {
            m_closestDistance = fDist;
            m_hasIntersection = true;
        }
    }
};
}

float Cortex::traceClosestHit(const BVH& bvhRef, const float3& origin, const float3& dir, float tMin, float tMax) const
{
    class RayHit : public IRay {
    public:
        float m_closest; bool m_found;
        RayHit(const float3& o, const float3& d, float tmin, float tmax)
        { m_vPos = o; m_vDir = normalize(d); m_fMin = tmin; m_fMax = tmax; m_closest = std::numeric_limits<float>::max(); m_found = false; }
        void notifyIntersection(float fDist, const ITraceableObject*, uint32_t) override
        { if (fDist >= m_fMin && fDist <= m_fMax && fDist < m_closest) { m_closest = fDist; m_found = true; } }
    } ray(origin, dir, tMin, tMax);
    bvhRef.trace(ray, 0);
    return ray.m_found ? ray.m_closest : 0.0f;
}

float3 Cortex::normalizedToWorld(const float3& normalizedPos)
{
    assert(m_pCortexBVH && "Cortex BVH must be initialized before normalizedToWorld");

    // Get bounding box center as ray origin
    const box3 boundingBox = m_pCortexBVH->getBox();
    const float3 center = boundingBox.center();

    if (length(normalizedPos) < 1e-6f) {
        return center;
    }

    const float3 direction = normalize(normalizedPos);
    const BVH& bvh = m_pCortexBVH->updateAndGetBVH();
    float hitDist = traceClosestHit(bvh, center, direction, 0.0f, std::numeric_limits<float>::max());
    if (hitDist <= 0.0f)
    {
        const float3 diagonal = boundingBox.diagonal();
        return center + normalizedPos * (diagonal * 0.5f);
    }

    float normalizedLength = length(normalizedPos);
    normalizedLength = std::min(1.0f, std::max(-1.0f, normalizedLength));
    return center + direction * (hitDist * normalizedLength);
}

float Cortex::distanceToCortex(const float3& originWorld, const float3& dirWorld)
{
    assert(m_pCortexBVH && "Cortex BVH must be initialized before distanceToCortex");

    const BVH& bvh = m_pCortexBVH->updateAndGetBVH();
    return traceClosestHit(bvh, originWorld, dirWorld, 0.0f, std::numeric_limits<float>::max());
}

float3 Cortex::baryToNormalized(uint32_t triangleIndex, const float3& barycentric) const
{
    if (!m_pTensionSphere)
        return float3(0,0,0);
    auto pMesh = m_pTensionSphere->getEdgeMesh();
    if (!pMesh)
        return float3(0,0,0);
    if (triangleIndex >= pMesh->getTriangleCount())
        return float3(0,0,0);

    const uint3 tri = pMesh->getTriangleVertices(triangleIndex);
    const float3 v0 = pMesh->getVertexPosition(tri.x);
    const float3 v1 = pMesh->getVertexPosition(tri.y);
    const float3 v2 = pMesh->getVertexPosition(tri.z);

    const float3 world = v0 * barycentric.x + v1 * barycentric.y + v2 * barycentric.z;

    // Find the direction of the ray from the center of the cube to the point of interest
    const box3 bbox = pMesh->getBox();
    const float3 center = bbox.center();
    const float3 half = (bbox.m_maxs - bbox.m_mins) * 0.5f;
    float3 dirWorld = world - center;

    // Convert to normalized cube space by scaling with half-extents
    const float eps = 1e-8f;
    float3 dN(
        (half.x > eps) ? (dirWorld.x / half.x) : 0.0f,
        (half.y > eps) ? (dirWorld.y / half.y) : 0.0f,
        (half.z > eps) ? (dirWorld.z / half.z) : 0.0f
    );

    // For the intersection with the unit cube surface, scale so the max abs component is 1
    float fMax = std::max(std::abs(dN.x), std::abs(dN.y));
    fMax = std::max(fMax, std::abs(dN.z));
    if (fMax < eps)
        return float3(0,0,0);
    return dN / fMax;
}