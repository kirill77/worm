#pragma once

#include <memory>
#include <vector>

class GPUMesh;

// this object is in 'gpu' and not in 'helpers' project because GPUWorld
// needs to know the declaration of this object
struct IVisObject
{
    virtual std::vector<std::shared_ptr<GPUMesh>> updateAndGetGpuMeshes() = 0;
}; 