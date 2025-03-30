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
    
    // Create cells based on faces
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
    return static_cast<uint32_t>(m_faces.size());
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

std::string TensionSphere::edgeKey(uint32_t v1, uint32_t v2) const
{
    // Ensure consistent ordering v1 < v2 for key generation
    if (v1 > v2) {
        std::swap(v1, v2);
    }
    
    // Create string key for edge
    return std::to_string(v1) + ":" + std::to_string(v2);
}

TensionSphere::Edge TensionSphere::createEdge(uint32_t startVertex, uint32_t endVertex)
{
    Edge edge;
    edge.startVertex = startVertex;
    edge.endVertex = endVertex;
    edge.leftFace = INVALID_INDEX;
    edge.rightFace = INVALID_INDEX;
    edge.startCW = INVALID_INDEX;
    edge.startCCW = INVALID_INDEX;
    edge.endCW = INVALID_INDEX;
    edge.endCCW = INVALID_INDEX;
    edge.leftCW = INVALID_INDEX;
    edge.leftCCW = INVALID_INDEX;
    edge.rightCW = INVALID_INDEX;
    edge.rightCCW = INVALID_INDEX;
    
    return edge;
}

uint32_t TensionSphere::addEdge(uint32_t startVertex, uint32_t endVertex)
{
    // Create new edge
    Edge edge = createEdge(startVertex, endVertex);
    uint32_t edgeIndex = static_cast<uint32_t>(m_edges.size());
    m_edges.push_back(edge);
    
    // Update vertex references to this edge
    if (m_vertices[startVertex].edgeIndex == INVALID_INDEX) {
        m_vertices[startVertex].edgeIndex = edgeIndex;
    }
    if (m_vertices[endVertex].edgeIndex == INVALID_INDEX) {
        m_vertices[endVertex].edgeIndex = edgeIndex;
    }
    
    // Add to edge map for faster lookups
    m_edgeMap[edgeKey(startVertex, endVertex)] = edgeIndex;
    
    return edgeIndex;
}

uint32_t TensionSphere::findEdge(uint32_t startVertex, uint32_t endVertex) const
{
    auto it = m_edgeMap.find(edgeKey(startVertex, endVertex));
    if (it != m_edgeMap.end()) {
        return it->second;
    }
    
    return INVALID_INDEX;
}

uint32_t TensionSphere::findOrCreateEdge(uint32_t startVertex, uint32_t endVertex)
{
    uint32_t edgeIndex = findEdge(startVertex, endVertex);
    if (edgeIndex == INVALID_INDEX) {
        edgeIndex = addEdge(startVertex, endVertex);
    }
    return edgeIndex;
}

void TensionSphere::linkEdges(uint32_t e1, uint32_t e2, bool startToStart, bool preserveWinding)
{
    Edge& edge1 = m_edges[e1];
    Edge& edge2 = m_edges[e2];
    
    // Determine which vertices to link
    uint32_t v1, v2;
    if (startToStart) {
        v1 = edge1.startVertex;
        v2 = edge2.startVertex;
    } else {
        v1 = edge1.startVertex;
        v2 = edge2.endVertex;
    }
    
    // Link edges based on winding order
    if (preserveWinding) {
        // Link e1 -> e2, same winding direction
        if (v1 == edge1.startVertex) {
            edge1.startCW = e2;
        } else {
            edge1.endCCW = e2;
        }
        
        if (v2 == edge2.startVertex) {
            edge2.startCCW = e1;
        } else {
            edge2.endCW = e1;
        }
    } else {
        // Link e1 -> e2, opposite winding direction
        if (v1 == edge1.startVertex) {
            edge1.startCCW = e2;
        } else {
            edge1.endCW = e2;
        }
        
        if (v2 == edge2.startVertex) {
            edge2.startCW = e1;
        } else {
            edge2.endCCW = e1;
        }
    }
}

void TensionSphere::attachEdgeToFace(uint32_t edgeIndex, uint32_t faceIndex, bool isLeftFace)
{
    // Get edge and update its face reference
    Edge& edge = m_edges[edgeIndex];
    if (isLeftFace) {
        edge.leftFace = faceIndex;
    } else {
        edge.rightFace = faceIndex;
    }
    
    // Update face's edge reference if not set
    Face& face = m_faces[faceIndex];
    if (face.edgeIndex == INVALID_INDEX) {
        face.edgeIndex = edgeIndex;
    }
}

void TensionSphere::completeEdgeLoop(uint32_t firstEdge, uint32_t lastEdge, uint32_t faceIndex)
{
    // Link edges within the face loop
    uint32_t e1 = firstEdge;
    uint32_t e2 = m_edges[e1].leftCW;
    
    while (e2 != INVALID_INDEX && e2 != lastEdge) {
        // Link e1 to e2 around the face
        m_edges[e1].leftCCW = e2;
        m_edges[e2].leftCW = e1;
        
        // Move to next edge
        e1 = e2;
        e2 = m_edges[e1].leftCW;
    }
    
    // Link last edge to complete the loop
    if (e2 == lastEdge) {
        m_edges[e1].leftCCW = e2;
        m_edges[e2].leftCW = e1;
        
        // Link lastEdge to firstEdge
        m_edges[lastEdge].leftCCW = firstEdge;
        m_edges[firstEdge].leftCW = lastEdge;
    }
}

uint32_t TensionSphere::addFace(uint32_t v1, uint32_t v2, uint32_t v3)
{
    // Create the face
    Face face;
    face.edgeIndex = INVALID_INDEX;
    face.area = 0.0;
    face.restArea = 0.0;
    
    uint32_t faceIndex = static_cast<uint32_t>(m_faces.size());
    m_faces.push_back(face);
    
    // Find or create edges
    uint32_t e1 = findOrCreateEdge(v1, v2);
    uint32_t e2 = findOrCreateEdge(v2, v3);
    uint32_t e3 = findOrCreateEdge(v3, v1);
    
    // Attach edges to face
    attachEdgeToFace(e1, faceIndex, true);
    attachEdgeToFace(e2, faceIndex, true);
    attachEdgeToFace(e3, faceIndex, true);
    
    // Link edges in cyclic order around face
    m_edges[e1].leftCCW = e2;
    m_edges[e2].leftCW = e1;
    
    m_edges[e2].leftCCW = e3;
    m_edges[e3].leftCW = e2;
    
    m_edges[e3].leftCCW = e1;
    m_edges[e1].leftCW = e3;
    
    // Set face's reference edge
    m_faces[faceIndex].edgeIndex = e1;
    
    // Link edges around vertices
    // For v1 (connects e3 -> e1)
    if (m_edges[e3].endVertex == v1) {
        m_edges[e3].endCCW = e1;
    } else {
        m_edges[e3].startCW = e1;
    }
    
    if (m_edges[e1].startVertex == v1) {
        m_edges[e1].startCW = e3;
    } else {
        m_edges[e1].endCCW = e3;
    }
    
    // For v2 (connects e1 -> e2)
    if (m_edges[e1].endVertex == v2) {
        m_edges[e1].endCCW = e2;
    } else {
        m_edges[e1].startCW = e2;
    }
    
    if (m_edges[e2].startVertex == v2) {
        m_edges[e2].startCW = e1;
    } else {
        m_edges[e2].endCCW = e1;
    }
    
    // For v3 (connects e2 -> e3)
    if (m_edges[e2].endVertex == v3) {
        m_edges[e2].endCCW = e3;
    } else {
        m_edges[e2].startCW = e3;
    }
    
    if (m_edges[e3].startVertex == v3) {
        m_edges[e3].startCW = e2;
    } else {
        m_edges[e3].endCCW = e2;
    }
    
    // Calculate and set area
    m_faces[faceIndex].area = calculateFaceArea(faceIndex);
    m_faces[faceIndex].restArea = m_faces[faceIndex].area;
    
    return faceIndex;
}

std::vector<uint32_t> TensionSphere::getFaceVertices(uint32_t faceIndex) const
{
    std::vector<uint32_t> vertices;
    
    if (faceIndex >= m_faces.size()) {
        return vertices;
    }
    
    const Face& face = m_faces[faceIndex];
    uint32_t firstEdge = face.edgeIndex;
    uint32_t currentEdge = firstEdge;
    
    do {
        const Edge& edge = m_edges[currentEdge];
        uint32_t v1 = edge.startVertex;
        
        vertices.push_back(v1);
        
        // Move to next edge around the face
        currentEdge = edge.leftCCW;
    } while (currentEdge != firstEdge && currentEdge != INVALID_INDEX && vertices.size() < 3);
    
    return vertices;
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
        double3 normal = calculateFaceNormal(i);
        
        // Force magnitude: positive tension pushes inward, negative pulls outward
        double forceMagnitude = tension * m_stiffness;
        
        // Apply force to each vertex (perpendicular to face)
        double3 forceVector = normal * forceMagnitude / 3.0;  // Divide by 3 to distribute to each vertex
        
        // Get vertices for this face and apply force
        std::vector<uint32_t> faceVertices = getFaceVertices(i);
        for (uint32_t vertexIndex : faceVertices) {
            m_vertices[vertexIndex].force += forceVector;
        }
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
        m_faces[i].area = calculateFaceArea(i);
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
        m_faces[i].area = calculateFaceArea(i);
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
    m_edges.clear();
    m_faces.clear();
    m_edgeMap.clear();
    
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
        m_vertices[i].edgeIndex = INVALID_INDEX;
    }
    
    // Define the 20 triangular faces of the icosahedron
    const uint32_t faceIndices[][3] = {
        {0, 8, 4}, {0, 4, 6}, {0, 6, 9}, {0, 9, 2}, {0, 2, 8},
        {1, 4, 8}, {1, 6, 4}, {1, 11, 6}, {1, 10, 11}, {1, 8, 10},
        {2, 7, 5}, {2, 9, 7}, {2, 5, 8}, {3, 5, 7}, {3, 7, 11},
        {3, 11, 10}, {3, 10, 5}, {4, 10, 8}, {5, 10, 8}, {6, 11, 9}
    };
    
    // Add faces with winged-edge data structure connections
    for (int i = 0; i < 20; ++i) {
        addFace(faceIndices[i][0], faceIndices[i][1], faceIndices[i][2]);
    }
}

void TensionSphere::subdivide(uint32_t level)
{
    for (uint32_t l = 0; l < level; ++l) {
        // Store the original faces and clear edge mapping
        std::vector<Face> originalFaces = m_faces;
        std::vector<Edge> originalEdges = m_edges;
        m_faces.clear();
        m_edgeMap.clear();
        
        // Keep track of midpoints for reuse
        std::unordered_map<std::string, uint32_t> midpointCache;
        
        // Process each original face
        for (uint32_t faceIdx = 0; faceIdx < originalFaces.size(); ++faceIdx) {
            // Get the three vertices of the face
            std::vector<uint32_t> vertices = getFaceVertices(faceIdx);
            if (vertices.size() != 3) continue;
            
            uint32_t v1 = vertices[0];
            uint32_t v2 = vertices[1];
            uint32_t v3 = vertices[2];
            
            // Get midpoints for each edge
            uint32_t mid1 = getMidpoint(v1, v2);
            uint32_t mid2 = getMidpoint(v2, v3);
            uint32_t mid3 = getMidpoint(v3, v1);
            
            // Create four new triangular faces
            addFace(v1, mid1, mid3);
            addFace(v2, mid2, mid1);
            addFace(v3, mid3, mid2);
            addFace(mid1, mid2, mid3);
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
    newVertex.edgeIndex = INVALID_INDEX;
    
    // Add the new vertex
    uint32_t newIndex = static_cast<uint32_t>(m_vertices.size());
    m_vertices.push_back(newVertex);
    
    // Store in cache and return
    midpointCache[key] = newIndex;
    return newIndex;
}

void TensionSphere::setupCellNeighbors()
{
    // Reset all neighbor lists
    for (auto& cell : m_cells) {
        cell.clearNeighbors();
    }
    
    // For each edge, if it connects two faces, they are neighbors
    for (const auto& edge : m_edges) {
        if (edge.leftFace != INVALID_INDEX && edge.rightFace != INVALID_INDEX) {
            // These faces share an edge, so they're neighbors
            m_cells[edge.leftFace].addNeighbor(edge.rightFace);
            m_cells[edge.rightFace].addNeighbor(edge.leftFace);
        }
    }
}

void TensionSphere::normalizeVertex(Vertex& v)
{
    v.position = normalize(v.position) * m_sphereRadius;
}

double TensionSphere::calculateFaceArea(uint32_t faceIndex) const
{
    // Get vertices of the face
    std::vector<uint32_t> vertices = getFaceVertices(faceIndex);
    if (vertices.size() != 3) return 0.0;
    
    const double3& p1 = m_vertices[vertices[0]].position;
    const double3& p2 = m_vertices[vertices[1]].position;
    const double3& p3 = m_vertices[vertices[2]].position;
    
    // Calculate vectors for two sides of the triangle
    double3 v1 = p2 - p1;
    double3 v2 = p3 - p1;
    
    // Cross product to get normal vector
    double3 crossProduct = cross(v1, v2);
    
    // Calculate area (half the length of the cross product vector)
    return 0.5 * length(crossProduct);
}

double3 TensionSphere::calculateFaceNormal(uint32_t faceIndex) const
{
    // Get vertices of the face
    std::vector<uint32_t> vertices = getFaceVertices(faceIndex);
    if (vertices.size() != 3) return double3(0.0, 0.0, 1.0);
    
    const double3& p1 = m_vertices[vertices[0]].position;
    const double3& p2 = m_vertices[vertices[1]].position;
    const double3& p3 = m_vertices[vertices[2]].position;
    
    // Calculate vectors for two sides of the triangle
    double3 v1 = p2 - p1;
    double3 v2 = p3 - p1;
    
    // Cross product to get normal vector
    double3 normal = cross(v1, v2);
    
    // Normalize
    return normalize(normal);
}
