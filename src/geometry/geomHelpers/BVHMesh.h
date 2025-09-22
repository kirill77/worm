#pragma once

#include "geometry/BVH/BVH.h"
#include "geometry/mesh/Mesh.h"

class BVHMesh : public ITraceableObject
{
public:
    BVHMesh(std::shared_ptr<Mesh> pMesh);

    // ITraceableObject interface
    virtual box3 getBox() const override;
    virtual box3 getSubObjectBox(uint32_t uSubObj) const override;
    virtual void trace(IRay& ray, uint32_t triangleIndex) const override;

    const BVH &updateAndGetBVH();
    
    std::shared_ptr<Mesh> getMesh() const { return m_pMesh; }

private:
    std::shared_ptr<Mesh> m_pMesh;
    BVH m_bvh;
    uint64_t m_cachedVersion = UINT64_MAX; // Invalid version initially
};

