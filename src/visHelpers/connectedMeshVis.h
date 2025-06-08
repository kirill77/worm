#pragma once

#include <memory>
#include <math/box.h>

class Window;
class ConnectedMesh;
class GPUMesh;
class GPUQueue;

struct ConnectedMeshVis
{
    ConnectedMeshVis(std::shared_ptr<Window>);

    void setConnectedBox(std::shared_ptr<box3> pBox) { m_pBox = pBox; }
    std::shared_ptr<box3> getConnectedBox() const { return m_pBox; }

    void setConnectedMesh(std::shared_ptr<ConnectedMesh> pMesh) { m_pMesh = pMesh; }

    std::shared_ptr<GPUMesh> getGPUMesh();
    void updateGPUMesh();

private:
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
    std::shared_ptr<box3> m_pBox;
};