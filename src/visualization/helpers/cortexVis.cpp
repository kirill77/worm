#include <memory>
#include "cortexVis.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/organelles/Cortex.h"
#include "geometry/mesh/edgeMesh.h"
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

    auto pMesh = m_pCortex->getMesh();

    // Get vertex count and triangle count
    uint32_t vertexCount = pMesh->getVertexCount();
    uint32_t triangleCount = pMesh->getTriangleCount();

    // Convert vertices to GPUMesh::Vertex format
    std::vector<GPUMesh::Vertex> gpuVertices;
    gpuVertices.reserve(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
    {
        GPUMesh::Vertex gpuVertex;
        convertVector(gpuVertex.vPos, pMesh->getVertexPosition(i));
        gpuVertices.push_back(gpuVertex);
    }

    // Convert triangles to GPU triangles
    std::vector<int3> gpuTriangles;
    gpuTriangles.reserve(triangleCount);
    for (uint32_t i = 0; i < triangleCount; ++i)
    {
        uint3 triangleVertices = pMesh->getTriangleVertices(i);
        gpuTriangles.push_back(int3(triangleVertices.x, triangleVertices.y, triangleVertices.z));
    }

    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
} 