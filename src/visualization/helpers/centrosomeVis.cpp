#include "CentrosomeVis.h"
#include "visualization/gpu/GPUQueue.h"
#include "visualization/gpu/GPUMeshNode.h"
#include "visualization/gpu/GPUMesh.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/Y_TuRC.h"
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
    m_pUnitCylinderGpuMesh = std::make_shared<GPUMesh>(pQueue->getDevice());
    
    // Create the static centrosome geometry once
    createCentrosomeGeometry();
    m_rootNode = GPUMeshNode(affine3::identity());
    // Initialize cached nodes once
    if (m_pUnitCylinderGpuMesh)
    {
        auto& children = m_rootNode.getChildren();
        children.clear();
        // child 0: X-axis node
        children.emplace_back(GPUMeshNode(affine3::identity()));
        children.back().addMesh(m_pUnitCylinderGpuMesh);
        // child 1: Y-axis node
        children.emplace_back(GPUMeshNode(affine3::identity()));
        children.back().addMesh(m_pUnitCylinderGpuMesh);
    }
    // nodes are initialized in constructor; nothing else to do here
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
    
    // Update transforms only
    m_rootNode.setTransform(centrosomeToWorld);

    updateCentriolesNodes();
    updateRingComplexNodes();

    return m_rootNode;
}

// Helper to construct a scaled cylinder orientation matrix (row-vector convention)
float3x3 CentrosomeVis::buildScaledCylinderMatrix(const float3& axisUnit, float length, float radius)
{
    float3 axis = normalize(axisUnit);
    if (length <= 0.0f) length = 0.0f;
    if (radius <= 0.0f) radius = 0.0f;

    float3 helper = (std::abs(axis.x) < 0.9f) ? float3(1, 0, 0) : float3(0, 1, 0);
    float3 radial1 = normalize(cross(axis, helper));
    float3 radial2 = cross(axis, radial1);

    return float3x3(
        radial1.x * radius, radial1.y * radius, radial1.z * radius,
        radial2.x * radius, radial2.y * radius, radial2.z * radius,
        axis.x * length,    axis.y * length,    axis.z * length
    );
}

void CentrosomeVis::updateCentriolesNodes()
{
    auto& children = m_rootNode.getChildren();
    if (children.size() < 2)
        return;

    // Centriole physical dimensions (micrometers) from EM literature; used for visualization.
    // Typical metazoan centriole: length ≈ 0.15 µm, radius ≈ 0.06 µm.
    const float fLengthMicroM = 0.15f; // µm
    const float fRadiusMicroM = 0.06f; // µm

    affine3 rotateToX = affine3::identity();
    rotateToX.m_linear = CentrosomeVis::buildScaledCylinderMatrix(float3(1,0,0), fLengthMicroM, fRadiusMicroM);
    children[0].setTransform(rotateToX);

    affine3 rotateToY = affine3::identity();
    rotateToY.m_linear = CentrosomeVis::buildScaledCylinderMatrix(float3(0,1,0), fLengthMicroM, fRadiusMicroM);
    children[1].setTransform(rotateToY);
}

// Visualize microtubules by scaling ring cylinders to match each MT's length
void CentrosomeVis::updateRingComplexNodes()
{
    auto& children = m_rootNode.getChildren();
    // Add/update additional cylinders for each ring complex (Y_TuRC)
    const auto& ringComplexes = m_pCentrosome->getRingComplexes();
    const size_t numRings = ringComplexes.size();

    // Ensure capacity: first two are centrioles
    const size_t needed = 2 + numRings;
    if (children.size() < needed && m_pUnitCylinderGpuMesh)
    {
        while (children.size() < needed)
        {
            GPUMeshNode n(affine3::identity());
            n.addMesh(m_pUnitCylinderGpuMesh);
            children.emplace_back(n);
        }
    }
    else if (children.size() > needed)
    {
        children.resize(needed);
    }

    // Ring cylinders are placed using Y_TuRC positions (in µm, relative to centrosome center)
    const float fDefaultStubMicroM = 0.04f;
    const float fRingCylRadiusMicroM = 0.01f;

    for (size_t i = 0; i < numRings; ++i)
    {
        auto& pRing = ringComplexes[i];
        assert(pRing && "Y_TuRC pointer should not be null");
        
        // For now, visualize as a single cylinder from origin to tip
        // TODO: Update to render as poly-line segments for bendable visualization
        float3 pos = pRing->getOrigin();
        float3 dir = pRing->getTipDirection();
        float3 tipPos = pRing->getTipPosition();

        affine3 ringXf = affine3::identity();
        float useLen = pRing->hasActiveMT() ? std::max(pRing->getMTLengthMicroM(), fDefaultStubMicroM) : fDefaultStubMicroM;
        ringXf.m_linear = CentrosomeVis::buildScaledCylinderMatrix(dir, useLen, fRingCylRadiusMicroM);
        // Shift by +0.5 * length along axis so the cylinder starts at the ring center and grows outward
        ringXf.m_translation = pos + dir * (0.5f * useLen);
        children[2 + i].setTransform(ringXf);
    }
}

void CentrosomeVis::createCentrosomeGeometry()
{
    // Create static cylinder geometry (single cylinder along Z-axis)
    // This will be reused for both X and Y orientations via transforms
    
    if (!m_pUnitCylinderGpuMesh)
    {
        assert(false);
        return;
    }

    const float radius = 1.0f;        // Cylinder radius
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
    m_pUnitCylinderGpuMesh->setGeometry(gpuVertices, gpuTriangles);
}
