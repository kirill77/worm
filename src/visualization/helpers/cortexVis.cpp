#include <memory>
#include "cortexVis.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/Cortex.h"
#include "geometry/edgeMesh/edgeMesh.h"
#include <stdexcept>

CortexVis::CortexVis(std::shared_ptr<Cortex> pCortex, GPUQueue *pQueue)
{
    m_pCortex = pCortex;
    m_pGPUMesh = std::make_shared<GPUMesh>(pQueue->getDevice());
}

std::shared_ptr<GPUMesh> CortexVis::updateAndGetGpuMesh()
{
    updateGPUMesh();
    return m_pGPUMesh;
}

void CortexVis::updateGPUMesh()
{
    if (!m_pCortex || !m_pGPUMesh)
    {
        return;
    }

    auto tensionSphere = m_pCortex->getTensionSphere();
    auto pMesh = tensionSphere.getEdgeMesh();

    // Get vertex count and face count
    uint32_t vertexCount = pMesh->getVertexCount();
    uint32_t faceCount = pMesh->getFaceCount();

    // Convert vertices to GPUMesh::Vertex format
    std::vector<GPUMesh::Vertex> gpuVertices;
    gpuVertices.reserve(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
    {
        GPUMesh::Vertex gpuVertex;
        convertVector(gpuVertex.vPos, pMesh->getVertexPosition(i));
        gpuVertices.push_back(gpuVertex);
    }

    // Convert faces to triangles
    std::vector<int3> gpuTriangles;
    gpuTriangles.reserve(faceCount);
    for (uint32_t i = 0; i < faceCount; ++i)
    {
        std::vector<uint32_t> faceVertices = pMesh->getFaceVertices(i);
        if (faceVertices.size() == 3)
        {
            gpuTriangles.push_back(int3(faceVertices[0], faceVertices[1], faceVertices[2]));
        }
    }

    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
} 