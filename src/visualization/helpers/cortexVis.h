#pragma once

#include "visualization/gpu/IVisObject.h"

class Organelle;
class GPUMesh;
struct GPUQueue;

struct CortexVis : public IVisObject
{
    CortexVis(std::shared_ptr<Organelle> pOrganelle, GPUQueue* pQueue);

    virtual GPUMeshNode updateAndGetMeshNode() override;

private:
    void updateGPUMesh();
    std::shared_ptr<Organelle> m_pOrganelle;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};