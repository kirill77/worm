#include "ForceGenerator.h"
#include "geometry/mesh/EdgeMesh.h"

EdgeSpringForce::EdgeSpringForce(IFaceBody& body, double springConstant)
    : m_body(body)
    , m_springConstant(springConstant)
{
    const uint32_t edgeCount = m_body.m_pMesh->getEdgeCount();
    m_edgeRestLengths.reserve(edgeCount);
    
    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto edge = m_body.m_pMesh->getEdge(e);
        float3 pos1 = m_body.m_pMesh->getVertexPosition(edge.first);
        float3 pos2 = m_body.m_pMesh->getVertexPosition(edge.second);
        double restLength = length(pos2 - pos1);
        m_edgeRestLengths.push_back(restLength);
    }
}

void EdgeSpringForce::apply(double dt)
{
    (void)dt;
    const uint32_t edgeCount = m_body.m_pMesh->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = m_body.m_pMesh->getEdge(e);
        double3 pa = double3(m_body.m_pMesh->getVertexPosition(ab.first));
        double3 pb = double3(m_body.m_pMesh->getVertexPosition(ab.second));
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

void EdgeDampingForce::apply(double dt)
{
    (void)dt;
    const uint32_t edgeCount = m_body.m_pMesh->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = m_body.m_pMesh->getEdge(e);
        double3 pa = double3(m_body.m_pMesh->getVertexPosition(ab.first));
        double3 pb = double3(m_body.m_pMesh->getVertexPosition(ab.second));
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


