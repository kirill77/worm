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
#include "VisObjectFactory.h"
#include "chemistry/ProteinWiki.h"
#include <windows.h>
#include "chemistry/StringDict.h"
#include "utils/log/ILog.h"
#include "visualization/gpu/DirectXHelpers.h"

bool VisEngine::initialize(std::shared_ptr<Organism> pOrganism)
{
    // Store the organism for later use
    m_pOrganism = pOrganism;

    // Initialize protein interaction data
    ProteinWiki::Initialize();

    // Create visualization window
    m_pWindow = std::make_shared<Window>();
    if (!m_pWindow->createWindowDevicAndSwapChain("Worm Simulation")) {
        LOG_ERROR("Failed to create window");
        return false;
    }

    // Create GPU world for visualization
    m_pGpuWorld = std::make_shared<GPUWorld>(m_pWindow, m_pWindow->getSwapChain()->getGPUQueue());

    // Initialize GPU text
    m_gpuText = std::make_unique<GPUText>(m_pGpuWorld->getFont());

    // Initialize camera UI
    m_cameraUI.attachToCamera(m_pGpuWorld->getCamera());

    // Initialize GPU stats
    m_gpuStats = std::make_unique<GPUStats>(m_pWindow->getDevice());

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

    // Check for space key press to toggle pause (ignore repeats)
    const UIState& uiState = m_pWindow->getCurrentUIState();
    const bool bIgoreRepeats = true;
    if (uiState.isPressed(VK_SPACE, bIgoreRepeats)) {
        m_bPaused = !m_bPaused;
    }

    if (!m_bPaused || // don't simulate if paused
        m_pWorld->getCurrentTime() < 5) { // need at least one step
        m_pWorld->simulateStep(fDtSec);
    }

    // Update GPU meshes
    updateGpuMeshes();

    // Update text with current simulation time and pause status
    if (m_bPaused) {
        m_gpuText->printf("%.2lf sec [PAUSED] - Press SPACE to resume", m_pWorld->getCurrentTime());
    } else {
        m_gpuText->printf("%.2lf sec - Press SPACE to pause", m_pWorld->getCurrentTime());
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
            m_cameraUI.fitWorldBoxToView();
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

void VisEngine::shutdown()
{
    if (m_pWindow && m_pWindow->getSwapChain()) {
        m_pWindow->getSwapChain()->getGPUQueue()->flush();
    }
}
