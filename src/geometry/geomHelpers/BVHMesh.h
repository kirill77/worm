#pragma once

#include "geometry/BVH/BVH.h"
#include "geometry/mesh/Mesh.h"
#include <cassert>

class BVHMesh : public ITraceableObject
{
public:
    BVHMesh(std::shared_ptr<Mesh> pMesh);

    // ITraceableObject interface
    virtual box3 getBox() const override;
    virtual box3 getSubObjectBox(uint32_t uSubObj) const override;
    virtual void trace(IRay& ray, uint32_t triangleIndex) const override;

    const BVH &getBVH() const {
        assert(m_pMesh && "BVHMesh must have a valid Mesh");
        assert(m_debugVersion == m_pMesh->getVersion() && "BVHMesh BVH is out of sync with Mesh version");
        return m_bvh;
    }

    // Rebuild BVH to match current Mesh topology
    void rebuildForCurrentMesh();
    
    std::shared_ptr<Mesh> getMesh() const { return m_pMesh; }

private:
    std::shared_ptr<Mesh> m_pMesh;
    BVH m_bvh;
    uint64_t m_debugVersion = 0;
};

