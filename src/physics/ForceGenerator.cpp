#include "ForceGenerator.h"
#include "geometry/mesh/TriangleMesh.h"
#include "geometry/mesh/Edges.h"

EdgeSpringForce::EdgeSpringForce(PhysicsMesh& body, double springConstant)
    : m_body(body)
    , m_springConstant(springConstant)
{
    auto pEdges = m_body.m_pMesh->getOrCreateEdges();
    const uint32_t edgeCount = pEdges->getEdgeCount();
    m_edgeRestLengths.reserve(edgeCount);
    
    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto edge = pEdges->getEdge(e);
        float3 pos1 = m_body.m_pMesh->getVertices()->getVertexPosition(edge.first);
        float3 pos2 = m_body.m_pMesh->getVertices()->getVertexPosition(edge.second);
        double restLength = length(pos2 - pos1);
        m_edgeRestLengths.push_back(restLength);
    }
}

void EdgeSpringForce::apply()
{
    auto pEdges = m_body.m_pMesh->getOrCreateEdges();
    const uint32_t edgeCount = pEdges->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = pEdges->getEdge(e);
        double3 pa = double3(m_body.m_pMesh->getVertices()->getVertexPosition(ab.first));
        double3 pb = double3(m_body.m_pMesh->getVertices()->getVertexPosition(ab.second));
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double L0 = m_edgeRestLengths[e];
        double3 f = -m_springConstant * (L - L0) * n;
        m_body.getVertex(ab.first).m_vForce += (-f);
        m_body.getVertex(ab.second).m_vForce += f;
    }
}

void EdgeDampingForce::apply()
{
    auto pEdges = m_body.m_pMesh->getOrCreateEdges();
    const uint32_t edgeCount = pEdges->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = pEdges->getEdge(e);
        double3 pa = double3(m_body.m_pMesh->getVertices()->getVertexPosition(ab.first));
        double3 pb = double3(m_body.m_pMesh->getVertices()->getVertexPosition(ab.second));
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double3 relV = m_body.getVertex(ab.second).m_vVelocity - m_body.getVertex(ab.first).m_vVelocity;
        double relAlong = dot(relV, n);
        double3 f = -m_dampingCoeff * relAlong * n;
        m_body.getVertex(ab.first).m_vForce += (-f);
        m_body.getVertex(ab.second).m_vForce += f;
    }
}


