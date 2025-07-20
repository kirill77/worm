#include "CentrosomeVis.h"
#include "visualization/gpu/GPUQueue.h"
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

std::shared_ptr<GPUMesh> CentrosomeVis::updateAndGetGpuMesh()
{
    updateGPUMesh();
    return m_pGPUMesh;
}

void CentrosomeVis::updateGPUMesh()
{
    if (!m_pCentrosome || !m_pGPUMesh)
    {
        assert(false);
        return;
    }

    // Get the cell from the organelle
    auto pCell = m_pCentrosome->getCell();
    if (!pCell)
    {
        assert(false);
        return;
    }
    // Get the cortex BVH for coordinate conversion
    auto pCortexBVH = pCell->getCortexBVH();
    if (!pCortexBVH)
    {
        assert(false);
        return;
    }

    // Convert normalized position to world coordinates
    const float3& normalizedPosition = m_pCentrosome->getNormalizedPosition();  // Position in normalized coordinates (-1, 1)
    float3 position = pCortexBVH->normalizedToWorld(normalizedPosition);  // Convert to world coordinates

    const float radius = 0.1f;        // Cylinder radius
    const float length = 0.8f;        // Cylinder length
    const int segments = 8;           // Number of segments around cylinder
    const float microtubuleLength = 2.0f;  // Length of microtubules
    const int microtubulesPerSegment = 2;  // Number of microtubules per cylinder segment
    
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
            float3 center = position + axis * (length * 0.5f * (end == 0 ? -1.0f : 1.0f));
            
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
    
    // Helper function to add microtubules sprouting from cylinder surface
    auto addMicrotubules = [&](const float3& axis, int startVertexOffset) {
        // Create two perpendicular vectors to the axis
        float3 perpVec1, perpVec2;
        if (abs(axis.x) < 0.9f) {
            perpVec1 = normalize(cross(axis, float3(1, 0, 0)));
        } else {
            perpVec1 = normalize(cross(axis, float3(0, 1, 0)));
        }
        perpVec2 = normalize(cross(axis, perpVec1));
        
        int currentVertexOffset = startVertexOffset;
        
        // Add microtubules along the cylinder length
        for (int lengthStep = 0; lengthStep < 3; ++lengthStep) {
            float t = (lengthStep / 2.0f) - 0.5f; // -0.5 to 0.5
            float3 cylinderCenter = position + axis * (length * t);
            
            // Add microtubules around the circumference
            for (int i = 0; i < segments; ++i) {
                for (int mt = 0; mt < microtubulesPerSegment; ++mt) {
                    float angle = 2.0f * 3.14159f * i / segments;
                    float angleOffset = (mt * 0.5f) / microtubulesPerSegment; // Small offset for multiple MTs
                    float finalAngle = angle + angleOffset;
                    
                    // Direction from cylinder surface outward
                    float3 radialDir = perpVec1 * cos(finalAngle) + perpVec2 * sin(finalAngle);
                    
                    // Add some randomness to make it more natural
                    float3 randomOffset = float3(
                        (lengthStep * 0.1f - 0.1f), 
                        (i * 0.05f - 0.2f), 
                        (mt * 0.08f - 0.04f)
                    );
                    float3 finalDir = normalize(radialDir + randomOffset * 0.3f);
                    
                    // Starting point on cylinder surface
                    float3 startPoint = cylinderCenter + radialDir * radius;
                    float3 endPoint = startPoint + finalDir * microtubuleLength;
                    
                    // Create microtubule as a thin parallelepiped (rectangular box)
                    float3 microtubuleDir = normalize(endPoint - startPoint);
                    float microtubuleWidth = 0.03f;  // Width of the microtubule
                    float microtubuleHeight = 0.02f; // Height of the microtubule
                    
                    // Create two perpendicular vectors to the microtubule direction
                    float3 widthVec, heightVec;
                    if (abs(microtubuleDir.z) < 0.9f) {
                        widthVec = normalize(cross(microtubuleDir, float3(0, 0, 1)));
                    } else {
                        widthVec = normalize(cross(microtubuleDir, float3(1, 0, 0)));
                    }
                    heightVec = normalize(cross(microtubuleDir, widthVec));
                    
                    // Scale the perpendicular vectors
                    widthVec *= microtubuleWidth * 0.5f;
                    heightVec *= microtubuleHeight * 0.5f;
                    
                    // Create 8 vertices for the parallelepiped
                    float3 corners[8];
                    corners[0] = startPoint - widthVec - heightVec;  // Bottom-left-start
                    corners[1] = startPoint + widthVec - heightVec;  // Bottom-right-start
                    corners[2] = startPoint + widthVec + heightVec;  // Top-right-start
                    corners[3] = startPoint - widthVec + heightVec;  // Top-left-start
                    corners[4] = endPoint - widthVec - heightVec;    // Bottom-left-end
                    corners[5] = endPoint + widthVec - heightVec;    // Bottom-right-end
                    corners[6] = endPoint + widthVec + heightVec;    // Top-right-end
                    corners[7] = endPoint - widthVec + heightVec;    // Top-left-end
                    
                    // Add vertices
                    int baseIdx = currentVertexOffset;
                    for (int v = 0; v < 8; ++v) {
                        GPUMesh::Vertex vertex;
                        convertVector(vertex.vPos, corners[v]);
                        gpuVertices.push_back(vertex);
                    }
                    
                    // Create triangles for the 6 faces of the parallelepiped
                    // Front face (start)
                    gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 1, baseIdx + 2));
                    gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 2, baseIdx + 3));
                    
                    // Back face (end)
                    gpuTriangles.push_back(int3(baseIdx + 4, baseIdx + 6, baseIdx + 5));
                    gpuTriangles.push_back(int3(baseIdx + 4, baseIdx + 7, baseIdx + 6));
                    
                    // Bottom face
                    gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 4, baseIdx + 1));
                    gpuTriangles.push_back(int3(baseIdx + 1, baseIdx + 4, baseIdx + 5));
                    
                    // Top face
                    gpuTriangles.push_back(int3(baseIdx + 2, baseIdx + 6, baseIdx + 3));
                    gpuTriangles.push_back(int3(baseIdx + 3, baseIdx + 6, baseIdx + 7));
                    
                    // Left face
                    gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 3, baseIdx + 4));
                    gpuTriangles.push_back(int3(baseIdx + 3, baseIdx + 7, baseIdx + 4));
                    
                    // Right face
                    gpuTriangles.push_back(int3(baseIdx + 1, baseIdx + 5, baseIdx + 2));
                    gpuTriangles.push_back(int3(baseIdx + 2, baseIdx + 5, baseIdx + 6));
                    
                    currentVertexOffset += 8;
                }
            }
        }
        
        return currentVertexOffset;
    };
    
    // Create first cylinder along X-axis
    int verticesAdded = createCylinder(float3(1, 0, 0), 0);
    
    // Create second cylinder along Y-axis
    verticesAdded = createCylinder(float3(0, 1, 0), verticesAdded);
    
#if 0
    // Add microtubules sprouting from both cylinders
    verticesAdded = addMicrotubules(float3(1, 0, 0), verticesAdded);
    verticesAdded = addMicrotubules(float3(0, 1, 0), verticesAdded);
    
    // Add some extra microtubules at the intersection (higher PCM density)
    int intersectionMicrotubules = 12;
    for (int i = 0; i < intersectionMicrotubules; ++i) {
        float theta = 2.0f * 3.14159f * i / intersectionMicrotubules;
        float phi = 3.14159f * (i % 3) / 3.0f; // Vary elevation
        
        float3 direction = float3(
            sin(phi) * cos(theta),
            sin(phi) * sin(theta),
            cos(phi)
        );
        
        float3 startPoint = position + direction * (radius * 1.2f);
        float3 endPoint = startPoint + direction * microtubuleLength;
        
        // Create intersection microtubule as a parallelepiped
        float3 microtubuleDir = normalize(endPoint - startPoint);
        float microtubuleWidth = 0.03f;
        float microtubuleHeight = 0.02f;
        
        // Create two perpendicular vectors to the microtubule direction
        float3 widthVec, heightVec;
        if (abs(microtubuleDir.z) < 0.9f) {
            widthVec = normalize(cross(microtubuleDir, float3(0, 0, 1)));
        } else {
            widthVec = normalize(cross(microtubuleDir, float3(1, 0, 0)));
        }
        heightVec = normalize(cross(microtubuleDir, widthVec));
        
        // Scale the perpendicular vectors
        widthVec *= microtubuleWidth * 0.5f;
        heightVec *= microtubuleHeight * 0.5f;
        
        // Create 8 vertices for the parallelepiped
        float3 corners[8];
        corners[0] = startPoint - widthVec - heightVec;
        corners[1] = startPoint + widthVec - heightVec;
        corners[2] = startPoint + widthVec + heightVec;
        corners[3] = startPoint - widthVec + heightVec;
        corners[4] = endPoint - widthVec - heightVec;
        corners[5] = endPoint + widthVec - heightVec;
        corners[6] = endPoint + widthVec + heightVec;
        corners[7] = endPoint - widthVec + heightVec;
        
        // Add vertices
        int baseIdx = verticesAdded;
        for (int v = 0; v < 8; ++v) {
            GPUMesh::Vertex vertex;
            convertVector(vertex.vPos, corners[v]);
            gpuVertices.push_back(vertex);
        }
        
        // Create triangles for the 6 faces
        // Front face
        gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 1, baseIdx + 2));
        gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 2, baseIdx + 3));
        
        // Back face
        gpuTriangles.push_back(int3(baseIdx + 4, baseIdx + 6, baseIdx + 5));
        gpuTriangles.push_back(int3(baseIdx + 4, baseIdx + 7, baseIdx + 6));
        
        // Bottom face
        gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 4, baseIdx + 1));
        gpuTriangles.push_back(int3(baseIdx + 1, baseIdx + 4, baseIdx + 5));
        
        // Top face
        gpuTriangles.push_back(int3(baseIdx + 2, baseIdx + 6, baseIdx + 3));
        gpuTriangles.push_back(int3(baseIdx + 3, baseIdx + 6, baseIdx + 7));
        
        // Left face
        gpuTriangles.push_back(int3(baseIdx + 0, baseIdx + 3, baseIdx + 4));
        gpuTriangles.push_back(int3(baseIdx + 3, baseIdx + 7, baseIdx + 4));
        
        // Right face
        gpuTriangles.push_back(int3(baseIdx + 1, baseIdx + 5, baseIdx + 2));
        gpuTriangles.push_back(int3(baseIdx + 2, baseIdx + 5, baseIdx + 6));
        
        verticesAdded += 8;
    }
#endif
    
    // Update the GPU mesh geometry
    m_pGPUMesh->setGeometry(gpuVertices, gpuTriangles);
}
