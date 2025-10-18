#include "ForceGenerator.h"
#include "geometry/mesh/EdgeMesh.h"

void EdgeSpringForce::apply(IFaceBody& body, double dt)
{
    (void)dt;
    const uint32_t edgeCount = body.m_pMesh->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = body.m_pMesh->getEdge(e);
        double3 pa = double3(body.m_pMesh->getVertexPosition(ab.first));
        double3 pb = double3(body.m_pMesh->getVertexPosition(ab.second));
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double L0 = m_edgeRestLengths[e];
        double3 f = -m_springConstant * (L - L0) * n;
        body.getVertex(ab.first).m_vForce += (-f);
        body.getVertex(ab.second).m_vForce += f;
    }
}

void EdgeDampingForce::apply(IFaceBody& body, double dt)
{
    (void)dt;
    const uint32_t edgeCount = body.m_pMesh->getEdgeCount();
    if (edgeCount == 0) return;

    for (uint32_t e = 0; e < edgeCount; ++e)
    {
        auto ab = body.m_pMesh->getEdge(e);
        double3 pa = double3(body.m_pMesh->getVertexPosition(ab.first));
        double3 pb = double3(body.m_pMesh->getVertexPosition(ab.second));
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double3 relV = body.getVertex(ab.second).m_vVelocity - body.getVertex(ab.first).m_vVelocity;
        double relAlong = dot(relV, n);
        double3 f = -m_dampingCoeff * relAlong * n;
        body.getVertex(ab.first).m_vForce += (-f);
        body.getVertex(ab.second).m_vForce += f;
    }
}


