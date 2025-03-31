#include "pch.h"
#include "WingedEdgeMesh.h"
#include <algorithm>
#include <cmath>

// Constants
constexpr double PHI = 1.61803398874989484820; // Golden ratio for icosahedron

// Edge constructor
WingedEdgeMesh::Edge::Edge()
    : startVertex(INVALID_INDEX)
    , endVertex(INVALID_INDEX)
    , leftFace(INVALID_INDEX)
    , rightFace(INVALID_INDEX)
    , startCW(INVALID_INDEX)
    , startCCW(INVALID_INDEX)
    , endCW(INVALID_INDEX)
    , endCCW(INVALID_INDEX)
    , leftCW(INVALID_INDEX)
    , leftCCW(INVALID_INDEX)
    , rightCW(INVALID_INDEX)
    , rightCCW(INVALID_INDEX)
{
}

// Default constructor - creates an empty mesh
WingedEdgeMesh::WingedEdgeMesh()
{
}

// Constructor to create an icosahedron with optional subdivision
WingedEdgeMesh::WingedEdgeMesh(double radius, uint32_t subdivisionLevel)
{
    createIcosahedron(radius);
    
    if (subdivisionLevel > 0) {
        subdivide(subdivisionLevel);
    }
}

// Clear the mesh
void WingedEdgeMesh::clear()
{
    m_vertices.clear();
    m_edges.clear();
    m_faces.clear();
    m_edgeMap.clear();
}

double3 WingedEdgeMesh::getVertexPosition(uint32_t index) const
{
    if (index >= m_vertices.size()) {
        return double3(0.0);
    }
    return m_vertices[index].position;
}

void WingedEdgeMesh::setVertexPosition(uint32_t index, const double3& position)
{
    if (index < m_vertices.size()) {
        m_vertices[index].position = position;
    }
}

// Add a vertex to the mesh
uint32_t WingedEdgeMesh::addVertex(const double3& position)
{
    Vertex vertex;
    vertex.position = position;
    vertex.edgeIndex = INVALID_INDEX;
    
    uint32_t index = static_cast<uint32_t>(m_vertices.size());
    m_vertices.push_back(vertex);
    return index;
}

// Generate a unique key for an edge
std::string WingedEdgeMesh::edgeKey(uint32_t v1, uint32_t v2) const
{
    // Ensure consistent ordering v1 < v2 for key generation
    if (v1 > v2) {
        std::swap(v1, v2);
    }
    
    // Create string key for edge
    return std::to_string(v1) + ":" + std::to_string(v2);
}

// Create an edge with default values
WingedEdgeMesh::Edge WingedEdgeMesh::createEdge(uint32_t startVertex, uint32_t endVertex)
{
    Edge edge;
    edge.startVertex = startVertex;
    edge.endVertex = endVertex;
    // Other fields are initialized by Edge constructor
    
    return edge;
}

// Add an edge to the mesh and return its index
uint32_t WingedEdgeMesh::addEdge(uint32_t startVertex, uint32_t endVertex)
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

// Find an edge by its vertex indices
uint32_t WingedEdgeMesh::findEdge(uint32_t startVertex, uint32_t endVertex) const
{
    auto it = m_edgeMap.find(edgeKey(startVertex, endVertex));
    if (it != m_edgeMap.end()) {
        return it->second;
    }
    
    return INVALID_INDEX;
}

// Find or create an edge
uint32_t WingedEdgeMesh::findOrCreateEdge(uint32_t startVertex, uint32_t endVertex)
{
    uint32_t edgeIndex = findEdge(startVertex, endVertex);
    if (edgeIndex == INVALID_INDEX) {
        edgeIndex = addEdge(startVertex, endVertex);
    }
    return edgeIndex;
}

// Link two edges around a vertex
void WingedEdgeMesh::linkEdges(uint32_t e1, uint32_t e2, bool startToStart, bool preserveWinding)
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

// Attach an edge to a face
void WingedEdgeMesh::attachEdgeToFace(uint32_t edgeIndex, uint32_t faceIndex, bool isLeftFace)
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

// Complete an edge loop around a face
void WingedEdgeMesh::completeEdgeLoop(uint32_t firstEdge, uint32_t lastEdge, uint32_t faceIndex)
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

// Add a triangular face to the mesh
uint32_t WingedEdgeMesh::addFace(uint32_t v1, uint32_t v2, uint32_t v3)
{
    // Create the face
    Face face;
    face.edgeIndex = INVALID_INDEX;
    
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
    
    return faceIndex;
}

// Get all vertices of a face
std::vector<uint32_t> WingedEdgeMesh::getFaceVertices(uint32_t faceIndex) const
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

// Get neighbors of a face (faces that share an edge)
std::vector<uint32_t> WingedEdgeMesh::getFaceNeighbors(uint32_t faceIndex) const
{
    std::vector<uint32_t> neighbors;
    
    if (faceIndex >= m_faces.size()) {
        return neighbors;
    }
    
    const Face& face = m_faces[faceIndex];
    uint32_t firstEdge = face.edgeIndex;
    uint32_t currentEdge = firstEdge;
    
    do {
        const Edge& edge = m_edges[currentEdge];
        // If this edge is shared with another face, add that face as a neighbor
        uint32_t neighborFace = (edge.leftFace == faceIndex) ? edge.rightFace : edge.leftFace;
        
        if (neighborFace != INVALID_INDEX) {
            // Check if this neighbor is already in the list
            if (std::find(neighbors.begin(), neighbors.end(), neighborFace) == neighbors.end()) {
                neighbors.push_back(neighborFace);
            }
        }
        
        // Move to next edge around the face
        currentEdge = edge.leftCCW;
    } while (currentEdge != firstEdge && currentEdge != INVALID_INDEX);
    
    return neighbors;
}

// Calculate face area
double WingedEdgeMesh::calculateFaceArea(uint32_t faceIndex) const
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

// Calculate face normal
double3 WingedEdgeMesh::calculateFaceNormal(uint32_t faceIndex) const
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

// Create an icosahedron (regular 20-faced polyhedron)
void WingedEdgeMesh::createIcosahedron(double radius)
{
    // Clear any existing mesh data
    clear();
    
    // Calculate constants for icosahedron construction
    const double norm = std::sqrt(1.0 + PHI * PHI);
    const double a = 1.0 / norm * radius;  // Shorter distance from origin
    const double b = PHI / norm * radius;  // Longer distance from origin
    
    // Define the 12 vertices of the icosahedron
    std::vector<double3> positions = {
        double3(0, a, b), double3(0, a, -b), double3(0, -a, b), double3(0, -a, -b),
        double3(a, b, 0), double3(a, -b, 0), double3(-a, b, 0), double3(-a, -b, 0),
        double3(b, 0, a), double3(-b, 0, a), double3(b, 0, -a), double3(-b, 0, -a)
    };
    
    // Add vertices
    m_vertices.reserve(positions.size());
    for (const auto& pos : positions) {
        addVertex(pos);
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

// Subdivide the mesh for smoother surfaces
void WingedEdgeMesh::subdivide(uint32_t levels)
{
    for (uint32_t l = 0; l < levels; ++l) {
        // Store the original faces and their vertex information
        std::vector<Face> originalFaces = m_faces;
        
        // Extract vertex information for each face before clearing faces
        std::vector<std::vector<uint32_t>> faceVertices;
        faceVertices.reserve(originalFaces.size());
        
        for (uint32_t faceIdx = 0; faceIdx < originalFaces.size(); ++faceIdx) {
            faceVertices.push_back(getFaceVertices(faceIdx));
        }
        
        // Now clear faces and edge mapping
        m_faces.clear();
        m_edgeMap.clear();
        
        // Keep track of midpoints for reuse
        std::unordered_map<std::string, uint32_t> midpointCache;
        
        // Calculate average radius of sphere
        double totalDistance = 0.0;
        for (const auto& vertex : m_vertices) {
            totalDistance += length(vertex.position);
        }
        double sphereRadius = totalDistance / m_vertices.size();
        
        // Process each original face
        for (uint32_t faceIdx = 0; faceIdx < originalFaces.size(); ++faceIdx) {
            const std::vector<uint32_t>& vertices = faceVertices[faceIdx];
            if (vertices.size() != 3) continue;
            
            uint32_t v1 = vertices[0];
            uint32_t v2 = vertices[1];
            uint32_t v3 = vertices[2];
            
            // Get midpoints for each edge
            uint32_t mid1 = getMidpoint(v1, v2, midpointCache, sphereRadius);
            uint32_t mid2 = getMidpoint(v2, v3, midpointCache, sphereRadius);
            uint32_t mid3 = getMidpoint(v3, v1, midpointCache, sphereRadius);
            
            // Create four new triangular faces
            addFace(v1, mid1, mid3);
            addFace(v2, mid2, mid1);
            addFace(v3, mid3, mid2);
            addFace(mid1, mid2, mid3);
        }
    }
}

// Get or create a midpoint between two vertices
uint32_t WingedEdgeMesh::getMidpoint(uint32_t v1, uint32_t v2, 
                                     std::unordered_map<std::string, uint32_t>& midpointCache,
                                     double sphereRadius)
{
    // Ensure v1 < v2 for consistent key generation
    if (v1 > v2) {
        std::swap(v1, v2);
    }
    
    // Create a key for the edge
    std::string key = std::to_string(v1) + ":" + std::to_string(v2);
    
    // Check if we already have this midpoint
    auto it = midpointCache.find(key);
    if (it != midpointCache.end()) {
        return it->second;
    }
    
    // Calculate the midpoint
    double3 midpoint = (m_vertices[v1].position + m_vertices[v2].position) * 0.5;
    
    // Normalize to keep on sphere
    midpoint = normalize(midpoint) * sphereRadius;
    
    // Add the new vertex
    uint32_t newIndex = addVertex(midpoint);
    
    // Store in cache and return
    midpointCache[key] = newIndex;
    return newIndex;
}
