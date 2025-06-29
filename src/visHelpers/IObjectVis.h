#pragma once

#include <memory>

class GPUMesh;

struct IObjectVis
{
    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() = 0;
};

