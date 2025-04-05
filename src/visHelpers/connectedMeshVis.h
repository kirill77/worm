#pragma once

#include <memory>

class Window;
class ConnectedMesh;
class GPUMesh;

struct ConnectedMeshVis
{
    ConnectedMeshVis(std::shared_ptr<Window>, std::shared_ptr<ConnectedMesh> pMesh);

    std::shared_ptr<GPUMesh> getGPUMesh();
    void updateGPUMesh();

private:
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};