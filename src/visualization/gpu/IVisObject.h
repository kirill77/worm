#pragma once

#include <memory>

class GPUMesh;

// this object is in 'gpu' and not in 'helpers' project because GPUWorld
// needs to know the declaration of this object
struct IVisObject
{
    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() = 0;
}; 