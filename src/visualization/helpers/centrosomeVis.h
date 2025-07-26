#pragma once

#include "visualization/gpu/IObjectVis.h"
#include <memory>

struct Centrosome;
class GPUMesh;
struct GPUQueue;

struct CentrosomeVis : public IObjectVis
{
    CentrosomeVis(std::shared_ptr<Centrosome> pCentrosome, GPUQueue* pQueue);

    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() override;

private:
    void updateGPUMesh();
    std::shared_ptr<Centrosome> m_pCentrosome;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};

