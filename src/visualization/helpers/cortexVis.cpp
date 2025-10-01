#include <memory>
#include "cortexVis.h"
#include "visualization/gpu/GPUMeshNode.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/organelles/Organelle.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "chemistry/molecules/StringDict.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "geometry/vectors/affine.h"
#include <stdexcept>

CortexVis::CortexVis(std::shared_ptr<Organelle> pOrganelle, GPUQueue *pQueue)
{
    m_pOrganelle = pOrganelle;
    m_pGPUMesh = std::make_shared<GPUMesh>(pQueue->getDevice());
}

GPUMeshNode CortexVis::updateAndGetMeshNode()
{
    updateGPUMesh();
    // Cortex uses identity transform since it defines the coordinate system
    GPUMeshNode node(affine3::identity());
    if (m_pGPUMesh)
    {
        node.addMesh(m_pGPUMesh);
    }
    return node;
}

void CortexVis::updateGPUMesh()
{
    if (!m_pOrganelle || !m_pGPUMesh)
    {
        assert(false);
        return;
    }

    auto pCell = m_pOrganelle->getCell();
    if (!pCell)
    {
        assert(false);
        return;
    }

    // Fetch BVH from Cortex organelle
    auto pCortex = std::dynamic_pointer_cast<Cortex>(pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    if (!pCortex)
    {
        assert(false);
        return;
    }
    auto pBVHMesh = pCortex->getBVHMesh();
    if (!pBVHMesh)
    {
        assert(false);
        return;
    }

    auto pMesh = pBVHMesh->getMesh();
    if (!pMesh)
    {
        assert(false);
        return;
    }

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