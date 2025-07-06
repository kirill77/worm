#include "pch.h"
#include "TensionSphere.h"
#include <cmath>
#include <algorithm>
#include <cassert>
#include <set>

constexpr double PI = 3.14159265358979323846;

// TensionSphere implementation
TensionSphere::TensionSphere(uint32_t subdivisionLevel)
{
    // Create the base mesh with icosahedron and subdivisions
    m_pMesh = std::make_shared<EdgeMesh>(1, subdivisionLevel);
    
    // Initialize physics simulation data
    initializePhysics();
}

void TensionSphere::initializePhysics()
{
    const uint32_t vertexCount = m_pMesh->getVertexCount();
    
    // Initialize velocities to zero
    m_vertexVelocities.resize(vertexCount, double3(0, 0, 0));
    
    // Build edge connectivity by examining all vertex pairs
    // Since we don't have direct edge access, we'll build edges from face information
    std::set<std::pair<uint32_t, uint32_t>> uniqueEdges;
    
    const uint32_t faceCount = m_pMesh->getFaceCount();
    for (uint32_t faceIdx = 0; faceIdx < faceCount; ++faceIdx)
    {
        std::vector<uint32_t> faceVertices = m_pMesh->getFaceVertices(faceIdx);
        if (faceVertices.size() == 3) // Triangle face
        {
            // Add the three edges of this triangle
            for (int i = 0; i < 3; ++i)
            {
                uint32_t v1 = faceVertices[i];
                uint32_t v2 = faceVertices[(i + 1) % 3];
                
                // Ensure consistent edge ordering (smaller index first)
                if (v1 > v2) std::swap(v1, v2);
                uniqueEdges.insert({v1, v2});
            }
        }
    }
    
    // Convert set to vector and compute rest lengths
    m_edgeConnectivity.clear();
    m_edgeRestLengths.clear();
    
    for (const auto& edge : uniqueEdges)
    {
        m_edgeConnectivity.push_back(edge);
        
        // Compute rest length as current distance between vertices
        double3 pos1 = m_pMesh->getVertexPosition(edge.first);
        double3 pos2 = m_pMesh->getVertexPosition(edge.second);
        double restLength = length(pos2 - pos1) * 1.5;
        m_edgeRestLengths.push_back(restLength);
    }
}

void TensionSphere::computeSpringForces(std::vector<double3>& forces, double dt)
{
    // Clear forces
    std::fill(forces.begin(), forces.end(), double3(0, 0, 0));
    
    // Compute spring forces for each edge
    for (size_t edgeIdx = 0; edgeIdx < m_edgeConnectivity.size(); ++edgeIdx)
    {
        uint32_t v1 = m_edgeConnectivity[edgeIdx].first;
        uint32_t v2 = m_edgeConnectivity[edgeIdx].second;
        
        double3 pos1 = m_pMesh->getVertexPosition(v1);
        double3 pos2 = m_pMesh->getVertexPosition(v2);
        double3 vel1 = m_vertexVelocities[v1];
        double3 vel2 = m_vertexVelocities[v2];
        
        // Edge vector and current length
        double3 edgeVector = pos2 - pos1;
        double currentLength = length(edgeVector);
        
        if (currentLength > 1e-10) // Avoid division by zero
        {
            double3 edgeDir = edgeVector / currentLength;
            double restLength = m_edgeRestLengths[edgeIdx];
            
            // Spring force: F = -k * (currentLength - restLength) * direction
            double springDisplacement = currentLength - restLength;
            double3 springForce = -m_fSpringC * springDisplacement * edgeDir;
            
            // Damping force: F = -c * relative_velocity_along_edge
            double3 relativeVel = vel2 - vel1;
            double relativeVelAlongEdge = dot(relativeVel, edgeDir);
            double3 dampingForce = -m_fDampingCoeff * relativeVelAlongEdge * edgeDir;
            
            // Total force on this edge
            double3 totalForce = springForce + dampingForce;
            
            // Apply equal and opposite forces to the two vertices
            forces[v1] -= totalForce;
            forces[v2] += totalForce;
        }
    }
}

void TensionSphere::integrateMotion(const std::vector<double3>& forces, double dt)
{
    const uint32_t vertexCount = m_pMesh->getVertexCount();
    
    // Simple Euler integration (could be improved with Verlet or RK4)
    for (uint32_t i = 0; i < vertexCount; ++i)
    {
        // Assume unit mass for simplicity (F = ma, so a = F when m = 1)
        double3 acceleration = forces[i];
        
        // Update velocity: v = v + a * dt
        m_vertexVelocities[i] += acceleration * dt;
        
        // Update position: x = x + v * dt  
        double3 currentPos = m_pMesh->getVertexPosition(i);
        double3 newPos = currentPos + m_vertexVelocities[i] * dt;
        
        m_pMesh->setVertexPosition(i, newPos);
    }
}

void TensionSphere::makeTimeStep(double fDtSec)
{
    const uint32_t vertexCount = m_pMesh->getVertexCount();
    std::vector<double3> forces(vertexCount);
    
    // Compute spring forces
    computeSpringForces(forces, fDtSec);
    
    // Integrate motion
    integrateMotion(forces, fDtSec);
}