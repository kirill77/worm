#pragma once

#include <memory>
#include <math/box.h>

class GPUMesh;

struct IObjectVis
{
    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() = 0;
    virtual std::shared_ptr<box3> getConnectedBox() const = 0;
};

