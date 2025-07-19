#include "BVHMesh.h"
#include <limits>
#include <algorithm>

BVHMesh::BVHMesh(std::shared_ptr<Mesh> pMesh)
    : m_pMesh(pMesh)
{
    if (m_pMesh)
    {
        m_nSubObjects = m_pMesh->getTriangleCount();
    }
    else
    {
        m_nSubObjects = 0;
    }
}

box3 BVHMesh::getBox()
{
    if (!m_pMesh || m_pMesh->getVertexCount() == 0)
    {
        return box3::empty();
    }

    // Get first vertex to initialize bounds
    double3 firstVertex = m_pMesh->getVertexPosition(0);
    float3 mins = toFloat3(firstVertex);
    float3 maxs = toFloat3(firstVertex);

    // Expand bounds to include all vertices
    for (uint32_t i = 1; i < m_pMesh->getVertexCount(); ++i)
    {
        float3 vertex = toFloat3(m_pMesh->getVertexPosition(i));
        mins = min(mins, vertex);
        maxs = max(maxs, vertex);
    }

    return box3(mins, maxs);
}

box3 BVHMesh::getSubObjectBox(uint32_t uSubObj)
{
    if (!m_pMesh || uSubObj >= m_pMesh->getTriangleCount())
    {
        return box3::empty();
    }

    // Get triangle vertices
    uint3 triangle = m_pMesh->getTriangleVertices(uSubObj);
    float3 v0 = toFloat3(m_pMesh->getVertexPosition(triangle.x));
    float3 v1 = toFloat3(m_pMesh->getVertexPosition(triangle.y));
    float3 v2 = toFloat3(m_pMesh->getVertexPosition(triangle.z));

    // Calculate bounding box of triangle
    float3 mins = min(min(v0, v1), v2);
    float3 maxs = max(max(v0, v1), v2);

    return box3(mins, maxs);
}

void BVHMesh::trace(IRay& ray)
{
    if (!m_pMesh)
        return;

    const float EPSILON = 1e-8f;

    // Test ray against each triangle using Möller-Trumbore algorithm
    for (uint32_t triangleIndex = 0; triangleIndex < m_pMesh->getTriangleCount(); ++triangleIndex)
    {
        // Get triangle vertices
        uint3 triangle = m_pMesh->getTriangleVertices(triangleIndex);
        float3 v0 = toFloat3(m_pMesh->getVertexPosition(triangle.x));
        float3 v1 = toFloat3(m_pMesh->getVertexPosition(triangle.y));
        float3 v2 = toFloat3(m_pMesh->getVertexPosition(triangle.z));

        // Möller-Trumbore ray-triangle intersection algorithm
        float3 edge1 = v1 - v0;
        float3 edge2 = v2 - v0;
        float3 h = cross(ray.m_vDir, edge2);
        float a = dot(edge1, h);

        // Check if ray is parallel to triangle
        if (a > -EPSILON && a < EPSILON)
            continue;

        float f = 1.0f / a;
        float3 s = ray.m_vPos - v0;
        float u = f * dot(s, h);

        // Check if intersection point is outside triangle
        if (u < 0.0f || u > 1.0f)
            continue;

        float3 q = cross(s, edge1);
        float v = f * dot(ray.m_vDir, q);

        // Check if intersection point is outside triangle
        if (v < 0.0f || u + v > 1.0f)
            continue;

        // Calculate t to find intersection point
        float t = f * dot(edge2, q);

        // Check if intersection is within ray bounds
        if (t > EPSILON && t >= ray.m_fMin && t <= ray.m_fMax)
        {
            // Valid intersection found
            ray.notifyIntersection(t, this, triangleIndex);
        }
    }
}
