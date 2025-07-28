#include "CentrosomeVis.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMeshNode.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/Cell.h"
#include "geometry/geomHelpers/BVHMesh.h"
#include "geometry/vectors/vector.h"
#include <memory>
#include <vector>
#include <cmath>

CentrosomeVis::CentrosomeVis(std::shared_ptr<Centrosome> pCentrosome, GPUQueue* pQueue)
{
    m_pCentrosome = pCentrosome;
    m_pGPUMesh = std::make_shared<GPUMesh>(pQueue->getDevice());
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
    
    auto pCortexBVH = pCell->getCortexBVH();
    if (!pCortexBVH)
    {
        return GPUMeshNode(affine3::identity());
    }

    // Get the centrosome's transform and convert to world coordinates
    const float3& normalizedPosition = m_pCentrosome->getNormalizedPosition();
    float3 position = pCortexBVH->normalizedToWorld(normalizedPosition);
    
    // Get the full transform matrix from centrosome space to cell space
    const affine3& centrosomeToCell = m_pCentrosome->getToParentTransform();
    
    // Convert the centrosome's transform to world space coordinates
    affine3 centrosomeToWorld = affine3::identity();
    centrosomeToWorld.m_translation = position;
    
    // Update the mesh geometry (without transform)
    updateGPUMesh();
    
    // Create node with the calculated transform
    GPUMeshNode node(centrosomeToWorld);
    if (m_pGPUMesh)
    {
        node.addMesh(m_pGPUMesh);
    }
    return node;
}

void CentrosomeVis::updateGPUMesh()
{
    if (!m_pCentrosome || !m_pGPUMesh)
    {
        assert(false);
        return;
    }

    const float radius = 0.1f;        // Cylinder radius
    const float length = 0.8f;        // Cylinder length
    const int segments = 8;           // Number of segments around cylinder
    
    std::vector<GPUMesh::Vertex> gpuVertices;
    std::vector<int3> gpuTriangles;
    
    // Helper function to create a cylinder
    auto createCylinder = [&](const float3& axis, int vertexOffset) {
        // Create two perpendicular vectors to the axis
        float3 perpVec1, perpVec2;
        if (abs(axis.x) < 0.9f) {
            perpVec1 = normalize(cross(axis, float3(1, 0, 0)));
        } else {
            perpVec1 = normalize(cross(axis, float3(0, 1, 0)));
        }
        perpVec2 = normalize(cross(axis, perpVec1));
        
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
            int centerIdx = vertexOffset + end * (segments + 1);
            
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
            int idx1 = vertexOffset + 1 + i;                    // First ring, current vertex
            int idx2 = vertexOffset + 1 + next;                 // First ring, next vertex
            int idx3 = vertexOffset + (segments + 1) + 1 + i;    // Second ring, current vertex
            int idx4 = vertexOffset + (segments + 1) + 1 + next; // Second ring, next vertex
            
            // Two triangles per quad
            gpuTriangles.push_back(int3(idx1, idx3, idx2));
            gpuTriangles.push_back(int3(idx2, idx3, idx4));
        }
        
        return 2 * (segments + 1); // Return number of vertices added
    };

    // Create first cylinder along X-axis
    int verticesAdded = createCylinder(float3(1, 0, 0), 0);
    
    // Create second cylinder along Y-axis
    verticesAdded = createCylinder(float3(0, 1, 0), verticesAdded);
    
    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
}
