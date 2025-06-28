#pragma once

#include "IObjectVis.h"

struct Window;
class ConnectedMesh;
class GPUMesh;
struct GPUQueue;

struct ConnectedMeshVis : public IObjectVis
{
    ConnectedMeshVis(std::shared_ptr<Window>);

    void setConnectedBox(std::shared_ptr<box3> pBox) { m_pBox = pBox; }
    virtual std::shared_ptr<box3> getConnectedBox() const override { return m_pBox; }

    void setConnectedMesh(std::shared_ptr<ConnectedMesh> pMesh) { m_pMesh = pMesh; }

    virtual std::shared_ptr<GPUMesh> updateAndGetGpuMesh() override;

private:
    void updateGPUMesh();
    std::shared_ptr<ConnectedMesh> m_pMesh;
    std::shared_ptr<GPUMesh> m_pGPUMesh;
    std::shared_ptr<box3> m_pBox;
};