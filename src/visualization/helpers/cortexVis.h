#pragma once

#include "IObjectVis.h"

class Organelle;
class GPUMesh;
struct GPUQueue;

struct CortexVis : public IObjectVis
{
    CortexVis(std::shared_ptr<Organelle> pOrganelle, GPUQueue* pQueue);

    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() override;

private:
    void updateGPUMesh();
    std::shared_ptr<Organelle> m_pOrganelle;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};