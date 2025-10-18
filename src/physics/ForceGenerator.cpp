#include "ForceGenerator.h"

void EdgeSpringForce::apply(IBody& body, double dt)
{
    (void)dt;
    auto& E = body.edges();
    auto& N = body.nodes();
    if (E.edgeCount() == 0) return;

    for (uint32_t e = 0; e < E.edgeCount(); ++e)
    {
        auto ab = E.edge(e);
        double3 pa = N.getPosition(ab.first);
        double3 pb = N.getPosition(ab.second);
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double L0 = E.restLength(e);
        double3 f = -m_springConstant * (L - L0) * n;
        N.addForce(ab.first, -f);
        N.addForce(ab.second,  f);
    }
}

void EdgeDampingForce::apply(IBody& body, double dt)
{
    (void)dt;
    auto& E = body.edges();
    auto& N = body.nodes();
    if (E.edgeCount() == 0) return;

    for (uint32_t e = 0; e < E.edgeCount(); ++e)
    {
        auto ab = E.edge(e);
        double3 pa = N.getPosition(ab.first);
        double3 pb = N.getPosition(ab.second);
        double3 edgeVec = pb - pa;
        double L = length(edgeVec);
        if (L <= 1e-10) continue;
        double3 n = edgeVec / L;
        double3 relV = N.getVelocity(ab.second) - N.getVelocity(ab.first);
        double relAlong = dot(relV, n);
        double3 f = -m_dampingCoeff * relAlong * n;
        N.addForce(ab.first, -f);
        N.addForce(ab.second,  f);
    }
}


