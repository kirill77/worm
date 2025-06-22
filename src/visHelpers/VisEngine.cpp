#include "VisEngine.h"
#include "simulation/World.h"
#include "simulation/Organism.h"
#include "simulation/Cell.h"
#include "visualization/Window.h"
#include "visualization/GPUWorld.h"
#include "visualization/GPUStats.h"
#include "visualization/GPUText.h"
#include "visHelpers/ConnectedMeshVis.h"
#include "visHelpers/CameraUI.h"
#include "molecules/ProteinWiki.h"
#include "log/ILog.h"
#include "visualization/DirectXHelpers.h"

static std::shared_ptr<ConnectedMeshVis> createCortexVis(
    std::shared_ptr<Organism> pOrganism,
    std::shared_ptr<Window> pWindow)
{
    // get a connected mesh that shows how the cortex looks like
    auto pCell = pOrganism->getCells()[0];
    auto pCortex = pCell->getCortex();
    auto pConnectedMesh = pCortex->getTensionSphere().getConnectedMesh();

    std::shared_ptr<ConnectedMeshVis> pCortexVis = std::make_shared<ConnectedMeshVis>(pWindow);
    pCortexVis->setConnectedMesh(pConnectedMesh);
    return pCortexVis;
}

bool VisEngine::initialize(std::shared_ptr<Organism> pOrganism)
{
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
    m_gpuText->printf("Hello World!");

    // Create cortex visualization
    m_pCortexVis = createCortexVis(pOrganism, m_pWindow);
    m_pGpuWorld->addMesh(m_pCortexVis->getGPUMesh());

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

    m_cameraUI.setWorldBox(*m_pCortexVis->getConnectedBox());
    m_cameraUI.notifyNewUIState(m_pWindow->getCurrentUIState());

    // Simulate one step
    m_pWorld->simulateStep(fDtSec);

    // Render the visualization
    m_pCortexVis->updateGPUMesh();

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

void VisEngine::shutdown()
{
    if (m_pWindow && m_pWindow->getSwapChain()) {
        m_pWindow->getSwapChain()->getGPUQueue()->flush();
    }
}
