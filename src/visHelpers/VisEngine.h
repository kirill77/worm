#pragma once

#include <memory>
#include "visHelpers/CameraUI.h"
#include "visualization/GPUStats.h"
#include "visualization/GPUText.h"

// Forward declarations
struct World;
struct Organism;
struct Window;
struct GPUWorld;
struct ConnectedMeshVis;
struct IObjectVis;
class GPUMesh;

struct VisEngine
{
    bool initialize(std::shared_ptr<Organism> pOrganism);
    bool update(float fDtSec);
    void shutdown();

private:
    std::shared_ptr<World> m_pWorld;
    std::shared_ptr<Window> m_pWindow;
    std::shared_ptr<GPUWorld> m_pGpuWorld;
    struct ObjectData
    {
        std::shared_ptr<IObjectVis> m_pObject; // the visualized object
        std::shared_ptr<GPUMesh> m_pGpuMesh; // that object converted to GPUMesh
    };
    std::vector<ObjectData> m_pObjects;
    std::unique_ptr<GPUText> m_gpuText;
    std::unique_ptr<GPUStats> m_gpuStats;
    CameraUI m_cameraUI;
};

