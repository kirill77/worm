#include "PhysicsMesh.h"
#include "geometry/mesh/EdgeMesh.h"

PhysicsMesh::PhysicsMesh(std::shared_ptr<EdgeMesh> mesh)
    : m_pMesh(std::move(mesh))
{
    const uint32_t vertexCount = m_pMesh->getVertices()->getVertexCount();
    m_nodeData.resize(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
    {
        m_nodeData[i].m_vVelocity = double3(0, 0, 0);
        m_nodeData[i].m_vForce = double3(0, 0, 0);
        m_nodeData[i].m_fMass = 1.0;
    }
}

