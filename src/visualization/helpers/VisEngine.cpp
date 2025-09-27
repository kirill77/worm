#include "VisEngine.h"
#include "biology/simulation/World.h"
#include "biology/simulation/Organism.h"
#include "biology/simulation/CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "visualization/gpu/Window.h"
#include "visualization/gpu/GPUWorld.h"
#include "visualization/gpu/GPUMesh.h"
#include "visualization/gpu/GPUStats.h"
#include "visualization/gpu/GPUText.h"
#include "visualization/helpers/CortexVis.h"
#include "visualization/helpers/CameraUI.h"
#include "visualization/helpers/CameraTransition.h"
#include "VisObjectFactory.h"
#include "chemistry/interactions/InteractionsWiki.h"
#include <windows.h>
#include "chemistry/molecules/StringDict.h"
#include "utils/log/ILog.h"
#include "visualization/gpu/DirectXHelpers.h"
#include "visualization/gpu/GPUCamera.h"
#include "biology/organelles/Centrosome.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include <cmath>
#include <algorithm>
#include <filesystem>
#include "utils/fileUtils/fileUtils.h"

VisEngine::~VisEngine() = default;

bool VisEngine::initialize(std::shared_ptr<Organism> pOrganism)
{
    // Initialize logging to a file in data/simOutDebug (Debug) or data/simOutRelease (Release)
    {
        std::filesystem::path dataPath;
        #ifdef _DEBUG
        const char* kSimOutFolder = "data/simOutDebug";
        #else
        const char* kSimOutFolder = "data/simOutRelease";
        #endif
        if (FileUtils::findTheFolder(kSimOutFolder, dataPath)) {
            std::filesystem::create_directories(dataPath);
            const std::string logPath = (dataPath / "sim.log").string();
            ILog::create(logPath);
        }
    }

    // Store the organism for later use
    m_pOrganism = pOrganism;

    // Initialize protein interaction data
    InteractionsWiki::Initialize();

    // Create visualization window
    m_pWindow = std::make_shared<Window>();
    if (!m_pWindow->createWindowDevicAndSwapChain("Worm Simulation")) {
        LOG_ERROR("Failed to create window");
        return false;
    }

    // Create GPU world for visualization
    m_pGpuWorld = std::make_shared<GPUWorld>(m_pWindow, m_pWindow->getSwapChain()->getGPUQueue());

    // Initialize GPU text
    m_gpuText = std::make_shared<GPUText>(m_pGpuWorld->getFont());
    
    // Create simulation time display line (lives forever - 0 seconds)
    m_pSimTimeLineText = m_gpuText->createLine();

    // Initialize camera UI
    m_cameraUI.attachToCamera(m_pGpuWorld->getCamera());

    // Initialize GPU stats
    m_gpuStats = std::make_shared<GPUStats>(m_pWindow->getDevice());

    // Create world
    m_pWorld = std::make_shared<World>(pOrganism);

    return true;
}

bool VisEngine::update(float fDtSec)
{
    // Process window messages
    m_pWindow->processMessages();
    if (m_pWindow->shouldExit())
        return false;

    m_cameraUI.notifyNewUIState(m_pWindow->getCurrentUIState());

    // Process UI input (keyboard, etc.)
    processUIMessages();

    // (UI already processed above)

    if (!m_bPaused || // don't simulate if paused
        m_pWorld->getCurrentTime() < 5) { // need at least one step
        m_pWorld->simulateStep(fDtSec);
    }

    // Update GPU meshes
    updateGpuMeshes();

    // Update active camera transition if any
    if (m_pCameraTransition)
    {
        bool bContinue = m_pCameraTransition->update(fDtSec);
        if (!bContinue)
        {
            m_pCameraTransition.reset();
        }
    }

    // Update text with current simulation time and pause status
    if (m_bPaused) {
        m_pSimTimeLineText->printf("%.2lf sec [PAUSED] - Press SPACE to resume", m_pWorld->getCurrentTime());
    } else {
        m_pSimTimeLineText->printf("%.2lf sec - Press SPACE to pause", m_pWorld->getCurrentTime());
    }

    auto* pSwapChain = m_pWindow->getSwapChain();
    auto pGpuQueue = pSwapChain->getGPUQueue();
    auto pCmdList = pGpuQueue->beginRecording();

    // Draw meshes into the command list and get the combined bounding box
    box3 combinedBoundingBox = m_pGpuWorld->render(pSwapChain, pCmdList.Get());

    // Set the world box with the combined bounding box of all meshes
    if (!combinedBoundingBox.isempty())
    {
        m_cameraUI.setWorldBox(combinedBoundingBox);
        double fVolume = combinedBoundingBox.computeVolume();
        // if the volume change is too large - re-fit the camera
        if (fVolume != 0 && abs(fVolume - m_fPrevFittedVolume) / fVolume > 0.5)
        {
            // Use camera's own helper to fit the new world bounding box
            if (auto cam = m_pGpuWorld->getCamera())
            {
                cam->fitBoxToView(combinedBoundingBox);
            }
        }
        m_fPrevFittedVolume = fVolume;
    }

    // Draw text on top
    m_gpuText->render(pSwapChain, m_pGpuWorld->getSharedRootSignature(), pCmdList.Get());

    pGpuQueue->execute(pCmdList);

    // Present the frame
    ThrowIfFailed(pSwapChain->getSwapChain()->Present(1, 0));

    return true;
}

void VisEngine::updateGpuMeshes()
{
    // List of organelles to visualize
    static const std::vector<StringDict::ID> organellesToVisualize = {
        StringDict::ID::ORGANELLE_CORTEX,
        StringDict::ID::ORGANELLE_CENTROSOME
    };

    // Process all organelles that have visualization contexts
    auto cells = m_pOrganism->getCellSims();
    for (auto& cell : cells)
    {
        // Loop through all organelles we want to visualize
        for (StringDict::ID organelleId : organellesToVisualize)
        {
            auto pOrganelle = cell->getCell()->getOrganelle(organelleId);
            if (pOrganelle)
            {
                // Initialize organelle visualization if needed
                auto pVisObject = pOrganelle->getVisObject();
                if (!pVisObject)
                {
                    pVisObject = VisObjectFactory::createForOrganelle(pOrganelle, organelleId,
                        m_pWindow->getSwapChain()->getGPUQueue());
                    if (pVisObject)
                    {
                        m_pGpuWorld->addObject(pVisObject);
                    }
                }
            }
        }
    }
}

void VisEngine::processUIMessages()
{
    const UIState& uiState = m_pWindow->getCurrentUIState();
    const bool bIgnoreRepeats = true;

    // SPACE toggles pause (ignore repeats)
    if (uiState.isPressed(VK_SPACE, bIgnoreRepeats))
    {
        m_bPaused = !m_bPaused;
    }

    // Tab key toggles camera transition to look at the cell's Centrosome
    if (uiState.isPressed(VK_TAB, bIgnoreRepeats))
    {
        auto pCurrentCamera = m_pGpuWorld->getCamera();
        m_pCameraTransition = m_camFocuser.goToNextFocus(m_pOrganism, pCurrentCamera);
        if (m_pCameraTransition)
        {
            m_cameraUI.setFocusBox(m_pCameraTransition->getFocusBox());
            
            // Create and display focus text for 5 seconds
            std::string focusName = m_camFocuser.getLastFocusedOrganelleName();
            auto pFocusLine = m_gpuText->createLine();
            pFocusLine->printf("Focusing on: %s", focusName.c_str());
            pFocusLine->setLifeTime(5);
        }
    }
}



void VisEngine::shutdown()
{
    if (m_pWindow && m_pWindow->getSwapChain()) {
        m_pWindow->getSwapChain()->getGPUQueue()->flush();
    }
}
