#pragma once

#include <memory>
#include "visualization/helpers/CameraUI.h"
#include "visualization/gpu/GPUStats.h"
#include "visualization/gpu/GPUText.h"
#include "chemistry/StringDict.h"

// Forward declarations
class World;
class Organism;
struct Window;
struct GPUWorld;
struct CortexVis;
class Organelle;

struct VisEngine
{
    bool initialize(std::shared_ptr<Organism> pOrganism);
    bool update(float fDtSec);
    void shutdown();
    
    // Access to the World object
    std::shared_ptr<World> getWorld() const { return m_pWorld; }

private:
    void updateGpuMeshes();
    std::shared_ptr<Organism> m_pOrganism;
    std::shared_ptr<World> m_pWorld;
    std::shared_ptr<Window> m_pWindow;
    std::shared_ptr<GPUWorld> m_pGpuWorld;
    std::unique_ptr<GPUText> m_gpuText;
    std::unique_ptr<GPUStats> m_gpuStats;

    CameraUI m_cameraUI;
    double m_fPrevFittedVolume = 0;
    bool m_bPaused = true;
};

