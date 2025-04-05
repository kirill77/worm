#pragma once

#include <memory>

class Window;
class ConnectedMesh;
class GPUMesh;
class GPUQueue;

struct ConnectedMeshVis
{
    ConnectedMeshVis(std::shared_ptr<Window>);

    void setConnectedMesh(std::shared_ptr<ConnectedMesh> pMesh);
    std::shared_ptr<GPUMesh> getGPUMesh();
    void updateGPUMesh(GPUQueue& gpuQueue);

private:
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};