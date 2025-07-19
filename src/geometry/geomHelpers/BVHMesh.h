#pragma once

#include "geometry/BVH/ITraceableObject.h"
#include "geometry/mesh/Mesh.h"

class BVHMesh : public ITraceableObject
{
public:
    BVHMesh(std::shared_ptr<Mesh> pMesh);

    // ITraceableObject interface
    virtual box3 getBox() override;
    virtual box3 getSubObjectBox(uint32_t uSubObj) override;
    virtual void trace(IRay& ray) override;

private:
    std::shared_ptr<Mesh> m_pMesh;
};

