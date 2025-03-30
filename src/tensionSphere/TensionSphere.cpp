#include "pch.h"
#include "TensionSphere.h"
#include <cmath>
#include <unordered_map>
#include <string>
#include <algorithm>
#include <cassert>

// Golden ratio used for icosahedron construction
constexpr double PHI = 1.61803398874989484820;
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
    // Create the base icosahedron
    createIcosahedron();
    
    // Subdivide if requested
    if (subdivisionLevel > 0) {
        subdivide(subdivisionLevel);
    }
    
    // Create cells based on faces and setup neighbors
    m_cells.resize(m_faces.size());
    for (uint32_t i = 0; i < m_faces.size(); ++i) {
        m_cells[i].setArea(m_faces[i].area);
        // Store rest area for force calculations
        m_faces[i].restArea = m_faces[i].area;
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
    assert(index < m_vertices.size());
    return m_vertices[index].position;
}

void TensionSphere::makeTimeStep(double fDtSec)
{
    // Vector to store new tension values
    std::vector<double> newTensions(m_cells.size());
    
    // Calculate new tensions based on area scalers, tension scalers, and neighbor influences
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        SphereCell& cell = m_cells[i];
        
        // Area influence: tension increases as area decreases below desired area
        double targetArea = m_faces[i].restArea * cell.getAreaScaler();
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
    for (auto& vertex : m_vertices) {
        vertex.force = double3(0.0);
    }
    
    // Calculate force for each face based on tension and area
    for (uint32_t i = 0; i < m_faces.size(); ++i) {
        const Face& face = m_faces[i];
        const SphereCell& cell = m_cells[i];
        double tension = cell.getTension();
        
        // Skip faces with negligible tension
        if (std::abs(tension) < 1e-6) continue;
        
        // Get face normal vector (normalized)
        double3 normal = calculateFaceNormal(face);
        
        // Force magnitude: positive tension pushes inward, negative pulls outward
        double forceMagnitude = tension * m_stiffness;
        
        // Apply force to each vertex (perpendicular to face)
        double3 forceVector = normal * forceMagnitude / 3.0;  // Divide by 3 to distribute to each vertex
        
        // Add to force accumulators
        m_vertices[face.v1].force += forceVector;
        m_vertices[face.v2].force += forceVector;
        m_vertices[face.v3].force += forceVector;
    }
}

void TensionSphere::integrateMotion(double fDtSec)
{
    // Simple semi-implicit Euler integration
    for (auto& vertex : m_vertices) {
        // Update velocity (with damping)
        vertex.velocity = vertex.velocity * (1.0 - m_damping) + vertex.force * fDtSec;
        
        // Update position
        vertex.position += vertex.velocity * fDtSec;
    }
}

void TensionSphere::enforceSphericalConstraint()
{
    // Ensure all vertices stay on the sphere surface
    for (auto& vertex : m_vertices) {
        // Calculate distance from origin
        double distance = length(vertex.position);
        
        // Skip if already at desired radius (within tolerance)
        if (std::abs(distance - m_sphereRadius) < 1e-6) continue;
        
        // Normalize and scale to sphere radius
        if (distance > 1e-10) {
            // Normalize and scale to sphere radius
            vertex.position = normalize(vertex.position) * m_sphereRadius;
            
            // Project velocity onto tangent plane
            double dotProduct = dot(vertex.velocity, vertex.position);
            dotProduct /= (m_sphereRadius * m_sphereRadius);
            
            vertex.velocity -= vertex.position * dotProduct;
        }
    }
}

void TensionSphere::updateCellAreas()
{
    // Recalculate areas for all faces after vertex movement
    for (uint32_t i = 0; i < m_faces.size(); ++i) {
        m_faces[i].area = calculateFaceArea(m_faces[i]);
        m_cells[i].setArea(m_faces[i].area);
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
    for (auto& vertex : m_vertices) {
        vertex.velocity = double3(0.0);
        vertex.force = double3(0.0);
    }
    
    // Reset all faces to their rest areas
    for (uint32_t i = 0; i < m_faces.size(); ++i) {
        m_faces[i].area = calculateFaceArea(m_faces[i]);
        m_faces[i].restArea = m_faces[i].area;
    }
    
    // Set all cells to balanced state
    for (uint32_t i = 0; i < m_cells.size(); ++i) {
        // Set area based on the corresponding face
        m_cells[i].setArea(m_faces[i].area);
        
        // Reset scalers to 1.0 (balanced)
        m_cells[i].setAreaScaler(1.0);
        m_cells[i].setTensionScaler(1.0);
        
        // Reset tension to 0.0 (relaxed state)
        m_cells[i].setTension(0.0);
    }
}

void TensionSphere::createIcosahedron()
{
    // Clear existing data
    m_vertices.clear();
    m_faces.clear();
    
    // Calculate constants for icosahedron construction
    const double norm = std::sqrt(1.0 + PHI * PHI);
    const double a = 1.0 / norm;  // Shorter distance from origin
    const double b = PHI / norm;  // Longer distance from origin
    
    // Define the 12 vertices of the icosahedron
    std::vector<double3> positions = {
        double3(0, a, b), double3(0, a, -b), double3(0, -a, b), double3(0, -a, -b),
        double3(a, b, 0), double3(a, -b, 0), double3(-a, b, 0), double3(-a, -b, 0),
        double3(b, 0, a), double3(-b, 0, a), double3(b, 0, -a), double3(-b, 0, -a)
    };
    
    // Add vertices with zero velocity and force
    m_vertices.resize(positions.size());
    for (size_t i = 0; i < positions.size(); ++i) {
        m_vertices[i].position = positions[i];
        m_vertices[i].velocity = double3(0.0);
        m_vertices[i].force = double3(0.0);
    }
    
    // Define the 20 triangular faces of the icosahedron
    const uint32_t faceIndices[][3] = {
        {0, 8, 4}, {0, 4, 6}, {0, 6, 9}, {0, 9, 2}, {0, 2, 8},
        {1, 4, 8}, {1, 6, 4}, {1, 11, 6}, {1, 10, 11}, {1, 8, 10},
        {2, 7, 5}, {2, 9, 7}, {2, 5, 8}, {3, 5, 7}, {3, 7, 11},
        {3, 11, 10}, {3, 10, 5}, {4, 10, 8}, {5, 10, 8}, {6, 11, 9}
    };
    
    // Add faces with area calculation
    for (int i = 0; i < 20; ++i) {
        Face face;
        face.v1 = faceIndices[i][0];
        face.v2 = faceIndices[i][1];
        face.v3 = faceIndices[i][2];
        face.area = 0.0;  // Will be calculated in a moment
        face.restArea = 0.0;
        m_faces.push_back(face);
        
        // Calculate area
        m_faces[i].area = calculateFaceArea(m_faces[i]);
        m_faces[i].restArea = m_faces[i].area;
    }
}

void TensionSphere::subdivide(uint32_t level)
{
    for (uint32_t l = 0; l < level; ++l) {
        // Store the original faces
        std::vector<Face> originalFaces = m_faces;
        m_faces.clear();
        
        // Midpoint cache to avoid duplicate vertices
        std::unordered_map<std::string, uint32_t> midpointCache;
        
        // Subdivide each face into 4 new triangles
        for (const auto& face : originalFaces) {
            // Get the three midpoints
            uint32_t mid1 = getMidpoint(face.v1, face.v2);
            uint32_t mid2 = getMidpoint(face.v2, face.v3);
            uint32_t mid3 = getMidpoint(face.v3, face.v1);
            
            // Create four new triangular faces
            Face newFace1 = {face.v1, mid1, mid3, 0.0, 0.0};
            Face newFace2 = {face.v2, mid2, mid1, 0.0, 0.0};
            Face newFace3 = {face.v3, mid3, mid2, 0.0, 0.0};
            Face newFace4 = {mid1, mid2, mid3, 0.0, 0.0};
            
            // Calculate areas for the new faces
            newFace1.area = calculateFaceArea(newFace1);
            newFace2.area = calculateFaceArea(newFace2);
            newFace3.area = calculateFaceArea(newFace3);
            newFace4.area = calculateFaceArea(newFace4);
            
            // Set rest areas
            newFace1.restArea = newFace1.area;
            newFace2.restArea = newFace2.area;
            newFace3.restArea = newFace3.area;
            newFace4.restArea = newFace4.area;
            
            // Add the new faces
            m_faces.push_back(newFace1);
            m_faces.push_back(newFace2);
            m_faces.push_back(newFace3);
            m_faces.push_back(newFace4);
        }
    }
}

uint32_t TensionSphere::getMidpoint(uint32_t v1, uint32_t v2)
{
    // Ensure v1 < v2 for consistent key generation
    if (v1 > v2) {
        std::swap(v1, v2);
    }
    
    // Create a key for the edge
    std::string key = std::to_string(v1) + ":" + std::to_string(v2);
    
    // Check if we already have this midpoint
    static std::unordered_map<std::string, uint32_t> midpointCache;
    auto it = midpointCache.find(key);
    if (it != midpointCache.end()) {
        return it->second;
    }
    
    // Calculate the midpoint
    double3 midpoint = (m_vertices[v1].position + m_vertices[v2].position) * 0.5;
    
    // Normalize to keep on unit sphere
    midpoint = normalize(midpoint) * m_sphereRadius;
    
    // Create new vertex
    Vertex newVertex;
    newVertex.position = midpoint;
    newVertex.velocity = double3(0.0);
    newVertex.force = double3(0.0);
    
    // Add the new vertex
    uint32_t newIndex = static_cast<uint32_t>(m_vertices.size());
    m_vertices.push_back(newVertex);
    
    // Store in cache and return
    midpointCache[key] = newIndex;
    return newIndex;
}

void TensionSphere::setupCellNeighbors()
{
    // Cells are neighbors if their corresponding faces share an edge (two vertices)
    for (uint32_t i = 0; i < m_faces.size(); ++i) {
        for (uint32_t j = i + 1; j < m_faces.size(); ++j) {
            // Count shared vertices
            int sharedVertices = 0;
            
            // Check each vertex combination
            if (m_faces[i].v1 == m_faces[j].v1 || m_faces[i].v1 == m_faces[j].v2 || m_faces[i].v1 == m_faces[j].v3) sharedVertices++;
            if (m_faces[i].v2 == m_faces[j].v1 || m_faces[i].v2 == m_faces[j].v2 || m_faces[i].v2 == m_faces[j].v3) sharedVertices++;
            if (m_faces[i].v3 == m_faces[j].v1 || m_faces[i].v3 == m_faces[j].v2 || m_faces[i].v3 == m_faces[j].v3) sharedVertices++;
            
            // If exactly 2 vertices are shared, the faces are neighbors
            if (sharedVertices == 2) {
                m_cells[i].addNeighbor(j);
                m_cells[j].addNeighbor(i);
            }
        }
    }
}

void TensionSphere::normalizeVertex(Vertex& v)
{
    v.position = normalize(v.position) * m_sphereRadius;
}

double TensionSphere::calculateFaceArea(const Face& face) const
{
    const double3& p1 = m_vertices[face.v1].position;
    const double3& p2 = m_vertices[face.v2].position;
    const double3& p3 = m_vertices[face.v3].position;
    
    // Calculate vectors for two sides of the triangle
    double3 v1 = p2 - p1;
    double3 v2 = p3 - p1;
    
    // Cross product to get normal vector
    double3 crossProduct = cross(v1, v2);
    
    // Calculate area (half the length of the cross product vector)
    return 0.5 * length(crossProduct);
}

double3 TensionSphere::calculateFaceNormal(const Face& face) const
{
    const double3& p1 = m_vertices[face.v1].position;
    const double3& p2 = m_vertices[face.v2].position;
    const double3& p3 = m_vertices[face.v3].position;
    
    // Calculate vectors for two sides of the triangle
    double3 v1 = p2 - p1;
    double3 v2 = p3 - p1;
    
    // Cross product to get normal vector
    double3 normal = cross(v1, v2);
    
    // Normalize
    return normalize(normal);
}
