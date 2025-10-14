#pragma once

#include <memory>
#include "visualization/helpers/CameraUI.h"
#include "visualization/helpers/CamFocuser.h"
#include "chemistry/molecules/StringDict.h"

// Forward declarations
class World;
class Organism;
struct Window;
struct GPUWorld;
struct CortexVis;
class Organelle;
class CameraTransition;
struct GPUStats;
struct GPUText;
struct Line;

struct VisEngine
{
    // Initialize logging output directory and log file; call before creating objects that log in constructors
    bool initLog();
    bool initialize(std::shared_ptr<Organism> pOrganism);
    bool update(float fDtSec);
    void shutdown();
    ~VisEngine();
    
    // Access to the World object
    std::shared_ptr<World> getWorld() const { return m_pWorld; }

private:
    void updateGpuMeshes();
    void processUIMessages();

    std::shared_ptr<Organism> m_pOrganism;
    std::shared_ptr<World> m_pWorld;
    std::shared_ptr<Window> m_pWindow;
    std::shared_ptr<GPUWorld> m_pGpuWorld;
    std::shared_ptr<GPUText> m_gpuText;
    std::shared_ptr<GPUStats> m_gpuStats;

    // Simulation time display line
    std::shared_ptr<Line> m_pSimTimeLineText;

    CameraUI m_cameraUI;
    CamFocuser m_camFocuser;
    double m_fPrevFittedVolume = 0;
    bool m_bPaused = true;

    // Smooth camera transition controller (e.g., to focus on Centrosome)
    std::shared_ptr<CameraTransition> m_pCameraTransition;
};

