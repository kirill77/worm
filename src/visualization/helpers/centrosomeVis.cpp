#include "CentrosomeVis.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMeshNode.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Cortex.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "geometry/vectors/vector.h"
#include "geometry/vectors/matrix.h"
#include <memory>
#include <vector>
#include <cmath>

CentrosomeVis::CentrosomeVis(std::shared_ptr<Centrosome> pCentrosome, GPUQueue* pQueue)
{
    m_pCentrosome = pCentrosome;
    m_pGPUMesh = std::make_shared<GPUMesh>(pQueue->getDevice());
    
    // Create the static centrosome geometry once
    createCentrosomeGeometry();
}

GPUMeshNode CentrosomeVis::updateAndGetMeshNode()
{
    if (!m_pCentrosome)
    {
        return GPUMeshNode(affine3::identity());
    }

    // Calculate the centrosome's transform
    auto pCell = m_pCentrosome->getCell();
    if (!pCell)
    {
        return GPUMeshNode(affine3::identity());
    }
    
    // Get the centrosome's transform and convert to world coordinates
    const float3& normalizedPosition = m_pCentrosome->getNormalizedPosition();
    // Use Cortex method for mapping normalized position to world
    auto pCortex = std::dynamic_pointer_cast<Cortex>(pCell->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
    float3 position = pCortex ? pCortex->normalizedToWorld(normalizedPosition) : float3(0,0,0);
    
    // Get the full transform matrix from centrosome space to cell space
    const affine3& centrosomeToCell = m_pCentrosome->getToParentTransform();
    
    // Convert the centrosome's transform to world space coordinates
    affine3 centrosomeToWorld = affine3::identity();
    centrosomeToWorld.m_translation = position;
    
    // Create root node with the centrosome's world position
    GPUMeshNode rootNode(centrosomeToWorld);
    
    if (m_pGPUMesh)
    {
        // Get scale factor based on PCM radius
        float scaleFactorBase = m_pCentrosome->getPCMRadius() * 0.8f;
        
        // Create first child node: cylinder along X-axis (rotate Z->X)
        // Rotation: 90 degrees around Y-axis to align Z with X
        affine3 rotateToX = affine3::identity();
        rotateToX.m_linear = float3x3(
            0, 0, scaleFactorBase,   // X = old Z * scale
            0, scaleFactorBase, 0,   // Y = old Y * scale
            -scaleFactorBase, 0, 0   // Z = -old X * scale
        );
        GPUMeshNode xAxisNode(rotateToX);
        xAxisNode.addMesh(m_pGPUMesh);
        
        // Create second child node: cylinder along Y-axis (rotate Z->Y)  
        // Rotation: -90 degrees around X-axis to align Z with Y
        affine3 rotateToY = affine3::identity();
        rotateToY.m_linear = float3x3(
            scaleFactorBase, 0, 0,   // X = old X * scale
            0, 0, scaleFactorBase,   // Y = old Z * scale
            0, -scaleFactorBase, 0   // Z = -old Y * scale
        );
        GPUMeshNode yAxisNode(rotateToY);
        yAxisNode.addMesh(m_pGPUMesh);
        
        // Add both child nodes to the root
        rootNode.addChild(std::move(xAxisNode));
        rootNode.addChild(std::move(yAxisNode));
    }
    
    return rootNode;
}

void CentrosomeVis::createCentrosomeGeometry()
{
    // Create static cylinder geometry (single cylinder along Z-axis)
    // This will be reused for both X and Y orientations via transforms
    
    if (!m_pGPUMesh)
    {
        assert(false);
        return;
    }

    const float radius = 0.1f;        // Cylinder radius
    const float length = 1.0f;        // Cylinder length
    const int segments = 8;           // Number of segments around cylinder
    
    std::vector<GPUMesh::Vertex> gpuVertices;
    std::vector<int3> gpuTriangles;
    
    // Create a single cylinder along Z-axis (0,0,1)
    float3 axis = float3(0, 0, 1);
    
    // Create two perpendicular vectors to the Z-axis
    float3 perpVec1 = float3(1, 0, 0);  // X-axis
    float3 perpVec2 = float3(0, 1, 0);  // Y-axis
    
    // Create vertices for cylinder ends
    for (int end = 0; end < 2; ++end) {
        float3 center = axis * (length * 0.5f * (end == 0 ? -1.0f : 1.0f));
        
        // Center vertex for the end cap
        GPUMesh::Vertex centerVertex;
        convertVector(centerVertex.vPos, center);
        gpuVertices.push_back(centerVertex);
        
        // Ring vertices for the end cap
        for (int i = 0; i < segments; ++i) {
            float angle = 2.0f * 3.14159f * i / segments;
            float3 offset = perpVec1 * (radius * cos(angle)) + perpVec2 * (radius * sin(angle));
            
            GPUMesh::Vertex vertex;
            convertVector(vertex.vPos, center + offset);
            gpuVertices.push_back(vertex);
        }
    }
    
    // Create triangles for end caps
    for (int end = 0; end < 2; ++end) {
        int centerIdx = end * (segments + 1);
        
        for (int i = 0; i < segments; ++i) {
            int next = (i + 1) % segments;
            int idx1 = centerIdx + 1 + i;
            int idx2 = centerIdx + 1 + next;
            
            if (end == 0) {
                // First end cap (reverse winding)
                gpuTriangles.push_back(int3(centerIdx, idx2, idx1));
            } else {
                // Second end cap (normal winding)
                gpuTriangles.push_back(int3(centerIdx, idx1, idx2));
            }
        }
    }
    
    // Create triangles for cylinder sides
    for (int i = 0; i < segments; ++i) {
        int next = (i + 1) % segments;
        
        // Indices for the ring vertices (skip center vertices)
        int idx1 = 1 + i;                    // First ring, current vertex
        int idx2 = 1 + next;                 // First ring, next vertex
        int idx3 = (segments + 1) + 1 + i;    // Second ring, current vertex
        int idx4 = (segments + 1) + 1 + next; // Second ring, next vertex
        
        // Two triangles per quad
        gpuTriangles.push_back(int3(idx1, idx2, idx3));
        gpuTriangles.push_back(int3(idx2, idx4, idx3));
    }
    
    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
}
