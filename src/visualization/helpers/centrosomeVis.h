#pragma once

#include "visualization/gpu/IVisObject.h"
#include <memory>

struct Centrosome;
class GPUMesh;
struct GPUQueue;

struct CentrosomeVis : public IVisObject
{
    CentrosomeVis(std::shared_ptr<Centrosome> pCentrosome, GPUQueue* pQueue);

    virtual GPUMeshNode updateAndGetMeshNode() override;

private:
    void createCentrosomeGeometry();
    std::shared_ptr<Centrosome> m_pCentrosome;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};

