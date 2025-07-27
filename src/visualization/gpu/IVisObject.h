#pragma once

#include <memory>

class GPUMesh;

struct IVisObject
{
    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() = 0;
}; 