#include "pch.h"
#include "TensionSphere.h"
#include <cmath>
#include <algorithm>
#include <cassert>

constexpr double PI = 3.14159265358979323846;

// SphereCell implementation
SphereCell::SphereCell()
    : m_fArea(0.0)
    , m_fAreaScaler(1.0)
    , m_fTension(0.0)
    , m_fTensionScaler(1.0)
{
}

void SphereCell::addNeighbor(uint32_t neighborIndex)
{
    // Check if this neighbor already exists
    if (std::find(m_neighbors.begin(), m_neighbors.end(), neighborIndex) == m_neighbors.end()) {
        m_neighbors.push_back(neighborIndex);
    }
}

// TensionSphere implementation
TensionSphere::TensionSphere(uint32_t subdivisionLevel)
    : m_simulationTime(0.0)
    , m_stiffness(10.0)        // Default stiffness
    , m_damping(0.8)           // Default damping
    , m_sphereRadius(1.0)      // Unit sphere
{
    // Create the base mesh with icosahedron and subdivisions
    m_mesh = ConnectedMesh(m_sphereRadius, subdivisionLevel);
    
    // Initialize vertex physics data
    m_vertexData.resize(m_mesh.getVertexCount());
    
    // Initialize face and cell data
    uint32_t faceCount = m_mesh.getFaceCount();
    m_faceData.resize(faceCount);
    m_cells.resize(faceCount);
    
    // Calculate initial areas
    for (uint32_t i = 0; i < faceCount; ++i) {
        m_faceData[i].area = m_mesh.calculateFaceArea(i);
        m_faceData[i].restArea = m_faceData[i].area;
        m_cells[i].setArea(m_faceData[i].area);
    }
    
    // Setup neighbor relationships
    setupCellNeighbors();
    
    // Initialize all cells to balanced state
    resetToBalancedState();
}

uint32_t TensionSphere::getCellCount() const
{
    return static_cast<uint32_t>(m_cells.size());
}

SphereCell& TensionSphere::getCell(uint32_t index)
{
    assert(index < m_cells.size());
    return m_cells[index];
}

const SphereCell& TensionSphere::getCell(uint32_t index) const
{
    assert(index < m_cells.size());
    return m_cells[index];
}

double3 TensionSphere::getVertexPosition(uint32_t index) const
{
    return m_mesh.getVertexPosition(index);
}

void TensionSphere::makeTimeStep(double fDtSec)
{
    // Vector to store new tension values
    std::vector<double> newTensions(m_cells.size());
    
    // Calculate new tensions based on area scalers, tension scalers, and neighbor influences
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        SphereCell& cell = m_cells[i];
        
        // Area influence: tension increases as area decreases below desired area
        double targetArea = m_faceData[i].restArea * cell.getAreaScaler();
        double currentArea = cell.getArea();
        double areaFactor = (targetArea > 0.0) ? (targetArea - currentArea) / targetArea : 0.0;
        
        // Base tension from area difference, scaled by tension scaler
        double baseTension = areaFactor * cell.getTensionScaler();
        
        // Apply damping to current tension
        double tensionDampingFactor = 0.8;
        double weightedTension = cell.getTension() * tensionDampingFactor + 
                                 baseTension * (1.0 - tensionDampingFactor);
        
        // Gather neighbor tensions for diffusion
        const auto& neighbors = cell.getNeighbors();
        double totalNeighborTension = 0.0;
        for (uint32_t neighborIdx : neighbors) {
            totalNeighborTension += m_cells[neighborIdx].getTension();
        }
        
        // Calculate diffusion factor (how much tension spreads to neighbors)
        double diffusionRate = 0.2 * fDtSec;
        double selfWeight = 1.0 - diffusionRate * neighbors.size();
        
        // Final tension is weighted average of cell's tension and neighbors' tensions
        newTensions[i] = selfWeight * weightedTension;
        if (!neighbors.empty()) {
            double avgNeighborTension = totalNeighborTension / neighbors.size();
            newTensions[i] += diffusionRate * neighbors.size() * avgNeighborTension;
        }
    }
    
    // Apply the new tensions
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        m_cells[i].setTension(newTensions[i]);
    }
    
    // Now handle vertex physics:
    
    // 1. Calculate forces on vertices based on tensions and areas
    calculateForces();
    
    // 2. Integrate motion using forces
    integrateMotion(fDtSec);
    
    // 3. Enforce spherical constraint (all vertices must remain on the sphere)
    enforceSphericalConstraint();
    
    // 4. Update cell areas after vertex movement
    updateCellAreas();
    
    // Update simulation time
    m_simulationTime += fDtSec;
}

void TensionSphere::calculateForces()
{
    // Reset all force accumulators
    for (auto& vertex : m_vertexData) {
        vertex.force = double3(0.0);
    }
    
    // Calculate force for each face based on tension and area
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        const Face& face = m_faceData[i];
        const SphereCell& cell = m_cells[i];
        double tension = cell.getTension();
        
        // Skip faces with negligible tension
        if (std::abs(tension) < 1e-6) continue;
        
        // Get face normal vector (normalized)
        double3 normal = calculateFaceNormal(i);
        
        // Force magnitude: positive tension pushes inward, negative pulls outward
        double forceMagnitude = tension * m_stiffness;
        
        // Apply force to each vertex (perpendicular to face)
        double3 forceVector = normal * forceMagnitude / 3.0;  // Divide by 3 to distribute to each vertex
        
        // Get vertices for this face and apply force
        std::vector<uint32_t> faceVertices = m_mesh.getFaceVertices(i);
        for (uint32_t vertexIndex : faceVertices) {
            m_vertexData[vertexIndex].force += forceVector;
        }
    }
}

void TensionSphere::integrateMotion(double fDtSec)
{
    // Simple semi-implicit Euler integration
    for (uint32_t i = 0; i < m_vertexData.size(); ++i) {
        Vertex& vertex = m_vertexData[i];
        
        // Update velocity (with damping)
        vertex.velocity = vertex.velocity * (1.0 - m_damping) + vertex.force * fDtSec;
        
        // Update position
        double3 newPosition = m_mesh.getVertexPosition(i) + vertex.velocity * fDtSec;
        m_mesh.setVertexPosition(i, newPosition);
    }
}

void TensionSphere::enforceSphericalConstraint()
{
    // Ensure all vertices stay on the sphere surface
    for (uint32_t i = 0; i < m_vertexData.size(); ++i) {
        Vertex& vertex = m_vertexData[i];
        double3 position = m_mesh.getVertexPosition(i);
        
        // Calculate distance from origin
        double distance = length(position);
        
        // Skip if already at desired radius (within tolerance)
        if (std::abs(distance - m_sphereRadius) < 1e-6) continue;
        
        // Normalize and scale to sphere radius
        if (distance > 1e-10) {
            // Normalize and scale to sphere radius
            position = normalize(position) * m_sphereRadius;
            m_mesh.setVertexPosition(i, position);
            
            // Project velocity onto tangent plane
            double dotProduct = dot(vertex.velocity, position);
            dotProduct /= (m_sphereRadius * m_sphereRadius);
            
            vertex.velocity -= position * dotProduct;
        }
    }
}

void TensionSphere::updateCellAreas()
{
    // Recalculate areas for all faces after vertex movement
    for (uint32_t i = 0; i < m_faceData.size(); ++i) {
        m_faceData[i].area = m_mesh.calculateFaceArea(i);
        m_cells[i].setArea(m_faceData[i].area);
    }
}

double TensionSphere::getTotalTensionEnergy() const
{
    double totalEnergy = 0.0;
    for (const auto& cell : m_cells) {
        // Energy proportional to square of tension times area
        totalEnergy += cell.getTension() * cell.getTension() * cell.getArea();
    }
    return totalEnergy;
}

void TensionSphere::resetToBalancedState()
{
    // Initialize all vertices with zero velocity and forces
    for (auto& vertex : m_vertexData) {
        vertex.velocity = double3(0.0);
        vertex.force = double3(0.0);
    }
    
    // Reset all faces to their rest areas
    for (uint32_t i = 0; i < m_faceData.size(); ++i) {
        m_faceData[i].area = m_mesh.calculateFaceArea(i);
        m_faceData[i].restArea = m_faceData[i].area;
    }
    
    // Set all cells to balanced state
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        // Set area based on the corresponding face
        m_cells[i].setArea(m_faceData[i].area);
        
        // Reset scalers to 1.0 (balanced)
        m_cells[i].setAreaScaler(1.0);
        m_cells[i].setTensionScaler(1.0);
        
        // Reset tension to 0.0 (relaxed state)
        m_cells[i].setTension(0.0);
    }
}

void TensionSphere::setupCellNeighbors()
{
    // Reset all neighbor lists
    for (auto& cell : m_cells) {
        cell.clearNeighbors();
    }
    
    // For each face, get its neighbors from the mesh
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        std::vector<uint32_t> neighbors = m_mesh.getFaceNeighbors(i);
        for (uint32_t neighborIdx : neighbors) {
            m_cells[i].addNeighbor(neighborIdx);
        }
    }
}

double3 TensionSphere::calculateFaceNormal(uint32_t faceIndex) const
{
    return m_mesh.calculateFaceNormal(faceIndex);
}
