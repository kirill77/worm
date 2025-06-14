#include "connectedMeshVis.h"
#include "visualization/GPUMesh.h"
#include "visualization/Window.h"
#include "connectedMesh/connectedMesh.h"
#include <stdexcept>

ConnectedMeshVis::ConnectedMeshVis(std::shared_ptr<Window> pWindow)
{
    m_pBox = std::make_shared<box3>(box<float, 3>::empty());
    m_pGPUMesh = std::make_shared<GPUMesh>(pWindow->getDevice());
}

std::shared_ptr<GPUMesh> ConnectedMeshVis::getGPUMesh()
{
    return m_pGPUMesh;
}

void ConnectedMeshVis::updateGPUMesh()
{
    if (!m_pMesh || !m_pGPUMesh)
    {
        return;
    }

    // Get vertex count and face count
    uint32_t vertexCount = m_pMesh->getVertexCount();
    uint32_t faceCount = m_pMesh->getFaceCount();

    // Convert vertices to GPUMesh::Vertex format
    std::vector<GPUMesh::Vertex> gpuVertices;
    gpuVertices.reserve(vertexCount);
    for (uint32_t i = 0; i < vertexCount; ++i)
    {
        GPUMesh::Vertex gpuVertex;
        convertVector(gpuVertex.vPos, m_pMesh->getVertexPosition(i));
        *m_pBox |= gpuVertex.vPos;
        gpuVertices.push_back(gpuVertex);
    }

    // Convert faces to triangles
    std::vector<int3> gpuTriangles;
    gpuTriangles.reserve(faceCount);
    for (uint32_t i = 0; i < faceCount; ++i)
    {
        std::vector<uint32_t> faceVertices = m_pMesh->getFaceVertices(i);
        if (faceVertices.size() == 3)
        {
            gpuTriangles.push_back(int3(faceVertices[0], faceVertices[1], faceVertices[2]));
        }
    }

    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
} 