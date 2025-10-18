#pragma once

#include <memory>
#include <vector>
#include "physics/BodyInterfaces.h"
#include "geometry/mesh/edgeMesh.h"

// Adapter that exposes EdgeMesh + per-vertex state as IFaceBody
class SoftBodyMeshAdapter : public IFaceBody
{
public:
    explicit SoftBodyMeshAdapter(std::shared_ptr<EdgeMesh> mesh)
    {
        m_pMesh = std::move(mesh);
        const uint32_t vertexCount = m_pMesh->getVertexCount();
        m_nodeData.resize(vertexCount);
        for (uint32_t i = 0; i < vertexCount; ++i)
        {
            m_nodeData[i].m_vVelocity = double3(0, 0, 0);
            m_nodeData[i].m_vForce = double3(0, 0, 0);
            m_nodeData[i].m_fMass = 1.0;
        }
    }

    const INodeView& getVertex(uint32_t index) const override
    {
        return m_nodeData[index];
    }

    INodeView& getVertex(uint32_t index) override
    {
        return m_nodeData[index];
    }

private:
    std::vector<INodeView> m_nodeData;
};


