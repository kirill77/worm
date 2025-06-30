#pragma once

#include <memory>
#include "visualization/helpers/CameraUI.h"
#include "visualization/helpers/VisObjectContext.h"
#include "visualization/gpu/GPUStats.h"
#include "visualization/gpu/GPUText.h"

// Forward declarations
struct World;
struct Organism;
struct Window;
struct GPUWorld;
struct ConnectedMeshVis;

struct VisEngine
{
    bool initialize(std::shared_ptr<Organism> pOrganism);
    bool update(float fDtSec);
    void shutdown();

private:
    void updateGpuMeshes();
    std::shared_ptr<Organism> m_pOrganism;
    std::shared_ptr<World> m_pWorld;
    std::shared_ptr<Window> m_pWindow;
    std::shared_ptr<GPUWorld> m_pGpuWorld;
    std::unique_ptr<GPUText> m_gpuText;
    std::unique_ptr<GPUStats> m_gpuStats;
    CameraUI m_cameraUI;
};

