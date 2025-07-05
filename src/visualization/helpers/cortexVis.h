#pragma once

#include "IObjectVis.h"

struct Cortex;
class GPUMesh;
struct GPUQueue;

struct CortexVis : public IObjectVis
{
    CortexVis(std::shared_ptr<Cortex> pCortex, GPUQueue* pQueue);

    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() override;

private:
    void updateGPUMesh();
    std::shared_ptr<Cortex> m_pCortex;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};