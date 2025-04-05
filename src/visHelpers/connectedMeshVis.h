#pragma once

struct ConnectedMeshVis
{
    ConnectedMeshVis(std::shared_ptr<ConnecteMesh> pMesh);

    std::shared_ptr<GPUMesh> getGPUMesh();
    void updateGPUMesh();

private:
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
};