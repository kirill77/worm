#pragma once

#include "IObjectVis.h"

struct Window;
class ConnectedMesh;
class GPUMesh;
struct GPUQueue;

struct CortexVis : public IObjectVis
{
    CortexVis(GPUQueue* pQueue);

    void setConnectedMesh(std::shared_ptr<ConnectedMesh> pMesh) { m_pMesh = pMesh; }

    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() override;

private:
    void updateGPUMesh();
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};