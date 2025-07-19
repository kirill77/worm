#include "BVHMesh.h"
#include <limits>
#include <algorithm>

// Custom ray implementation for finding mesh intersections
class BVHMeshRay : public IRay
{
public:
    float m_closestDistance;
    bool m_hasIntersection;
    
    BVHMeshRay(const float3& pos, const float3& dir, float maxDist = std::numeric_limits<float>::max())
        : m_closestDistance(std::numeric_limits<float>::max())
        , m_hasIntersection(false)
    {
        m_vPos = pos;
        m_vDir = normalize(dir);
        m_fMin = 0.0f;
        m_fMax = maxDist;
    }
    
    virtual void notifyIntersection(float fDist, const ITraceableObject* pObject, uint32_t uSubObj) override
    {
        if (fDist >= m_fMin && fDist <= m_fMax && fDist < m_closestDistance)
        {
            m_closestDistance = fDist;
            m_hasIntersection = true;
        }
    }
};

BVHMesh::BVHMesh(std::shared_ptr<Mesh> pMesh)
    : m_pMesh(pMesh)
{
    m_nSubObjects = m_pMesh->getTriangleCount();
}

const BVH &BVHMesh::updateAndGetBVH()
{
    // Check if the mesh has been modified since last update
    uint64_t currentVersion = m_pMesh->getVersion();
    if (m_cachedVersion != currentVersion)
    {
        // Mesh has changed or BVH doesn't exist, create/rebuild BVH
        m_nSubObjects = m_pMesh->getTriangleCount();
        
        // Create new BVH and add this mesh to it
        m_bvh.accessObjects().clear();
        m_bvh.accessObjects().push_back(shared_from_this());
        m_bvh.rebuildHierarchy();
        
        m_cachedVersion = currentVersion;
    }
    
    return m_bvh;
}

box3 BVHMesh::getBox() const
{
    return m_pMesh->getBox();
}

box3 BVHMesh::getSubObjectBox(uint32_t uSubObj) const
{
    // Get triangle vertices
    uint3 triangle = m_pMesh->getTriangleVertices(uSubObj);
    float3 v0 = m_pMesh->getVertexPosition(triangle.x);
    float3 v1 = m_pMesh->getVertexPosition(triangle.y);
    float3 v2 = m_pMesh->getVertexPosition(triangle.z);

    // Calculate bounding box of triangle
    float3 mins = min(min(v0, v1), v2);
    float3 maxs = max(max(v0, v1), v2);

    return box3(mins, maxs);
}

void BVHMesh::trace(IRay& ray, uint32_t triangleIndex) const
{
    assert(m_cachedVersion == m_pMesh->getVersion());
    const float EPSILON = 1e-8f;

    // Get triangle vertices
    uint3 triangle = m_pMesh->getTriangleVertices(triangleIndex);
    float3 v0 = m_pMesh->getVertexPosition(triangle.x);
    float3 v1 = m_pMesh->getVertexPosition(triangle.y);
    float3 v2 = m_pMesh->getVertexPosition(triangle.z);

    // Möller-Trumbore ray-triangle intersection algorithm
    float3 edge1 = v1 - v0;
    float3 edge2 = v2 - v0;
    float3 h = cross(ray.m_vDir, edge2);
    float a = dot(edge1, h);

    // Check if ray is parallel to triangle
    if (a > -EPSILON && a < EPSILON)
        return;

    float f = 1.0f / a;
    float3 s = ray.m_vPos - v0;
    float u = f * dot(s, h);

    // Check if intersection point is outside triangle
    if (u < 0.0f || u > 1.0f)
        return;

    float3 q = cross(s, edge1);
    float v = f * dot(ray.m_vDir, q);

    // Check if intersection point is outside triangle
    if (v < 0.0f || u + v > 1.0f)
        return;

    // Calculate t to find intersection point
    float t = f * dot(edge2, q);

    // Check if intersection is within ray bounds
    if (t > EPSILON && t >= ray.m_fMin && t <= ray.m_fMax)
    {
        // Valid intersection found
        ray.notifyIntersection(t, this, triangleIndex);
    }
}

float3 BVHMesh::normalizedToWorld(const float3& normalizedPos)
{
    // Get the bounding box of the mesh
    box3 boundingBox = getBox();
    
    // Calculate the center of the bounding box
    float3 center = boundingBox.center();
    
    // Create direction vector from center towards the normalized position
    // If normalizedPos is at origin, we can't determine direction, so return center
    if (length(normalizedPos) < 1e-6f)
    {
        return center;
    }
    
    float3 direction = normalize(normalizedPos);
    
    // Create ray from center in the direction of the normalized point
    BVHMeshRay ray(center, direction);
    
    // Update BVH and trace the ray to find intersection with mesh
    const BVH& bvh = updateAndGetBVH();
    bvh.trace(ray, 0);
    
    if (!ray.m_hasIntersection)
    {
        // If no intersection found, fallback to simple mapping using bounding box
        float3 diagonal = boundingBox.diagonal();
        float3 worldPos = center + normalizedPos * (diagonal * 0.5f);
        return worldPos;
    }
    
    // Linear interpolation based on intersection distance
    // normalizedPos is in [-1,1] range, intersection is at boundary (length 1)
    float normalizedLength = length(normalizedPos);
    
    // Clamp to [-1,1] range
    normalizedLength = std::min(1.0f, std::max(-1.0f, normalizedLength));
    
    // Map from normalized space to world space
    float3 worldPos = center + direction * (ray.m_closestDistance * normalizedLength);
    
    return worldPos;
}
