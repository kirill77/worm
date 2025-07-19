#pragma once

#include "geometry/BVH/BVH.h"
#include "geometry/mesh/Mesh.h"

class BVHMesh : public ITraceableObject
{
public:
    BVHMesh(std::shared_ptr<Mesh> pMesh);

    // ITraceableObject interface
    virtual box3 getBox() override;
    virtual box3 getSubObjectBox(uint32_t uSubObj) override;
    virtual void trace(IRay& ray, uint32_t triangleIndex) override;

    const BVH &updateAndGetBVH();

private:
    std::shared_ptr<Mesh> m_pMesh;
    BVH m_bvh;
    uint64_t m_cachedVersion = UINT64_MAX; // Invalid version initially
};

