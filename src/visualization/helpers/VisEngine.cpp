#include "VisEngine.h"
#include "biology/World.h"
#include "biology/Organism.h"
#include "biology/Cell.h"
#include "biology/Cortex.h"
#include "visualization/gpu/Window.h"
#include "visualization/gpu/GPUWorld.h"
#include "visualization/gpu/GPUStats.h"
#include "visualization/gpu/GPUText.h"
#include "visualization/helpers/CortexVis.h"
#include "visualization/helpers/CameraUI.h"
#include "chemistry/ProteinWiki.h"
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

    // Simulate one step
    m_pWorld->simulateStep(fDtSec);

    // Update GPU meshes
    updateGpuMeshes();

    // Update text with current simulation time
    m_gpuText->printf("%.2lf sec", m_pWorld->getCurrentTime());

    auto* pSwapChain = m_pWindow->getSwapChain();
    auto pGpuQueue = pSwapChain->getGPUQueue();
    auto pCmdList = pGpuQueue->beginRecording();

    // Draw meshes into the command list
    m_pGpuWorld->render(pSwapChain, pCmdList.Get());

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

    // Initialize combined bounding box as empty
    box3 combinedBoundingBox = box3::empty();

    // Process all organelles that have visualization contexts
    auto cells = m_pOrganism->getCells();
    for (auto& cell : cells)
    {
        // Loop through all organelles we want to visualize
        for (StringDict::ID organelleId : organellesToVisualize)
        {
            auto pOrganelle = cell->getOrganelle(organelleId);
            if (pOrganelle)
            {
                // Initialize organelle visualization if needed
                auto pVisContext = pOrganelle->getVisObjectContext();
                if (!pVisContext)
                {
                    VisObjectContext::createForOrganelle(pOrganelle, organelleId,
                        m_pWindow->getSwapChain()->getGPUQueue());
                    pVisContext = pOrganelle->getVisObjectContext();
                }

                if (pVisContext && pVisContext->m_pObject)
                {
                    auto pGpuMesh = pVisContext->m_pObject->updateAndGetGpuMesh();
           
                    // Add mesh to the world if needed
                    if (pVisContext->m_pGpuMesh != pGpuMesh)
                    {
                        pVisContext->m_pGpuMesh = pGpuMesh;
                        m_pGpuWorld->addMesh(pVisContext->m_pGpuMesh);
                    }
                    
                    // Combine this mesh's bounding box with the overall bounding box
                    if (pGpuMesh)
                    {
                        const box3& meshBoundingBox = pGpuMesh->getBoundingBox();
                        combinedBoundingBox = combinedBoundingBox | meshBoundingBox;
                    }
                }
            }
        }
    }
    
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
}

void VisEngine::shutdown()
{
    if (m_pWindow && m_pWindow->getSwapChain()) {
        m_pWindow->getSwapChain()->getGPUQueue()->flush();
    }
}
