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
    
    // Validate the mesh structure
    validateMesh();
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
    // This function is no longer needed in the classic Winged-Edge structure
    // since we don't maintain the vertex-edge connectivity directly.
    // This is a no-op to avoid breaking existing code.
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
    
    // Attach edges to face (check if already attached to another face)
    Edge& edge1 = m_edges[e1];
    Edge& edge2 = m_edges[e2];
    Edge& edge3 = m_edges[e3];
    
    // For each edge, decide whether to use leftFace or rightFace
    // If leftFace is already assigned, use rightFace
    if (edge1.leftFace == INVALID_INDEX) {
        attachEdgeToFace(e1, faceIndex, true);
    } else {
        attachEdgeToFace(e1, faceIndex, false);
    }
    
    if (edge2.leftFace == INVALID_INDEX) {
        attachEdgeToFace(e2, faceIndex, true);
    } else {
        attachEdgeToFace(e2, faceIndex, false);
    }
    
    if (edge3.leftFace == INVALID_INDEX) {
        attachEdgeToFace(e3, faceIndex, true);
    } else {
        attachEdgeToFace(e3, faceIndex, false);
    }
    
    // Link edges in cyclic order around face
    // Determine if this is a leftFace or rightFace for each edge
    bool e1IsLeft = (edge1.leftFace == faceIndex);
    bool e2IsLeft = (edge2.leftFace == faceIndex);
    bool e3IsLeft = (edge3.leftFace == faceIndex);
    
    // Link in proper CCW order for left faces, CW order for right faces
    if (e1IsLeft) {
        edge1.leftCCW = e2;
        if (e2IsLeft) {
            edge2.leftCW = e1;
        } else {
            edge2.rightCCW = e1;
        }
    } else {
        edge1.rightCW = e2;
        if (e2IsLeft) {
            edge2.leftCCW = e1;
        } else {
            edge2.rightCW = e1;
        }
    }
    
    if (e2IsLeft) {
        edge2.leftCCW = e3;
        if (e3IsLeft) {
            edge3.leftCW = e2;
        } else {
            edge3.rightCCW = e2;
        }
    } else {
        edge2.rightCW = e3;
        if (e3IsLeft) {
            edge3.leftCCW = e2;
        } else {
            edge3.rightCW = e2;
        }
    }
    
    if (e3IsLeft) {
        edge3.leftCCW = e1;
        if (e1IsLeft) {
            edge1.leftCW = e3;
        } else {
            edge1.rightCCW = e3;
        }
    } else {
        edge3.rightCW = e1;
        if (e1IsLeft) {
            edge1.leftCCW = e3;
        } else {
            edge1.rightCW = e3;
        }
    }
    
    // Set face's reference edge
    m_faces[faceIndex].edgeIndex = e1;
    
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
    
    // Return empty list if this face doesn't have a valid edge
    if (firstEdge == INVALID_INDEX || firstEdge >= m_edges.size()) {
        return vertices;
    }
    
    // Track visited edges to prevent infinite loops
    std::vector<uint32_t> visitedEdges;
    uint32_t currentEdge = firstEdge;
    
    // Check if this face is a left or right face for this edge
    bool isLeftFace = (m_edges[currentEdge].leftFace == faceIndex);
    
    // Safety counter to prevent infinite loops in case of inconsistent mesh
    int safetyCounter = 0;
    const int MAX_EDGES = 10; // Triangular faces have 3 edges, but be safe
    
    do {
        // Add current edge to visited list
        visitedEdges.push_back(currentEdge);
        
        const Edge& edge = m_edges[currentEdge];
        
        // For a left face, we traverse CCW and add start vertex
        // For a right face, we traverse CW and add end vertex
        uint32_t nextEdge;
        
        if (isLeftFace) {
            vertices.push_back(edge.startVertex);
            nextEdge = edge.leftCCW;
        } else {
            vertices.push_back(edge.endVertex);
            nextEdge = edge.rightCW;
        }
        
        // Check for invalid next edge
        if (nextEdge == INVALID_INDEX) {
            break;
        }
        
        // Check if we've already visited this edge (cycle detection)
        if (std::find(visitedEdges.begin(), visitedEdges.end(), nextEdge) != visitedEdges.end()) {
            break;
        }
        
        currentEdge = nextEdge;
        
        // Safety check
        safetyCounter++;
        if (safetyCounter >= MAX_EDGES) {
            break;
        }
        
    } while (currentEdge != firstEdge && vertices.size() < 3);
    
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
    
    // Return empty list if this face doesn't have a valid edge
    if (firstEdge == INVALID_INDEX || firstEdge >= m_edges.size()) {
        return neighbors;
    }
    
    // Track visited edges to prevent infinite loops
    std::vector<uint32_t> visitedEdges;
    uint32_t currentEdge = firstEdge;
    
    // Check if this face is a left or right face for this edge
    bool isLeftFace = (m_edges[currentEdge].leftFace == faceIndex);
    
    // Safety counter to prevent infinite loops in case of inconsistent mesh
    int safetyCounter = 0;
    const int MAX_EDGES = 10; // Faces shouldn't have more than this many edges
    
    do {
        // Add current edge to visited list
        visitedEdges.push_back(currentEdge);
        
        const Edge& edge = m_edges[currentEdge];
        
        // If this edge is shared with another face, add that face as a neighbor
        uint32_t neighborFace;
        uint32_t nextEdge;
        
        if (isLeftFace) {
            neighborFace = edge.rightFace;
            nextEdge = edge.leftCCW;
        } else {
            neighborFace = edge.leftFace;
            nextEdge = edge.rightCW;
        }
        
        if (neighborFace != INVALID_INDEX) {
            // Check if this neighbor is already in the list
            if (std::find(neighbors.begin(), neighbors.end(), neighborFace) == neighbors.end()) {
                neighbors.push_back(neighborFace);
            }
        }
        
        // Check for invalid next edge
        if (nextEdge == INVALID_INDEX) {
            break;
        }
        
        // Check if we've already visited this edge (cycle detection)
        if (std::find(visitedEdges.begin(), visitedEdges.end(), nextEdge) != visitedEdges.end()) {
            break;
        }
        
        currentEdge = nextEdge;
        
        // Safety check
        safetyCounter++;
        if (safetyCounter >= MAX_EDGES) {
            break;
        }
        
    } while (currentEdge != firstEdge);
    
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
    // Get vertices of the face in correct order
    std::vector<uint32_t> vertices = getFaceVertices(faceIndex);
    if (vertices.size() != 3) return double3(0.0, 0.0, 1.0);
    
    const double3& p1 = m_vertices[vertices[0]].position;
    const double3& p2 = m_vertices[vertices[1]].position;
    const double3& p3 = m_vertices[vertices[2]].position;
    
    // Calculate vectors for two sides of the triangle
    double3 v1 = p2 - p1;
    double3 v2 = p3 - p1;
    
    // Cross product to get normal vector (order matters for orientation)
    double3 normal = cross(v1, v2);
    
    // Normalize
    double len = length(normal);
    if (len > 1e-10) {
        normal = normal / len;
    } else {
        normal = double3(0.0, 0.0, 1.0);
    }
    
    return normal;
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
    
    // Validate the mesh structure after creation
    validateMesh();
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
        
        // Validate the mesh structure after each subdivision level
        validateMesh();
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
    
    // Check if we already have this midpoint in the cache
    auto it = midpointCache.find(key);
    if (it != midpointCache.end()) {
        return it->second;
    }
    
    // Check if an edge exists between these vertices
    uint32_t edgeIndex = findEdge(v1, v2);
    if (edgeIndex != INVALID_INDEX) {
        // If the edge exists, we need to keep track of its old data
        // to update it properly after the subdivision
        // Currently, we just create a new midpoint
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

// Validate mesh structure and fix inconsistencies
void WingedEdgeMesh::validateMesh()
{
    // Report statistics
    size_t vertexCount = m_vertices.size();
    size_t edgeCount = m_edges.size();
    size_t faceCount = m_faces.size();
    
    // 1. Verify each edge has both a left and right face (for a closed mesh)
    size_t edgesWithOneFace = 0;
    size_t edgesWithNoFace = 0;
    
    for (size_t i = 0; i < edgeCount; ++i) {
        Edge& edge = m_edges[i];
        
        // Check that edge vertices are valid indices
        if (edge.startVertex >= vertexCount || edge.endVertex >= vertexCount) {
            // Edge has invalid vertex indices - critical error
            continue;
        }
        
        // Check face connectivity
        bool hasLeftFace = (edge.leftFace != INVALID_INDEX && edge.leftFace < faceCount);
        bool hasRightFace = (edge.rightFace != INVALID_INDEX && edge.rightFace < faceCount);
        
        if (!hasLeftFace && !hasRightFace) {
            edgesWithNoFace++;
        } else if (hasLeftFace && !hasRightFace) {
            edgesWithOneFace++;
            
            // For a closed manifold mesh like an icosahedron, every edge should 
            // have exactly two adjacent faces. If only one is present, it might 
            // indicate an issue in face creation.
        } else if (!hasLeftFace && hasRightFace) {
            edgesWithOneFace++;
            
            // Swap left and right to ensure left is always populated
            std::swap(edge.leftFace, edge.rightFace);
            std::swap(edge.leftCW, edge.rightCW);
            std::swap(edge.leftCCW, edge.rightCCW);
        }
    }
    
    // 2. Verify face loop connectivity for each face
    for (size_t i = 0; i < faceCount; ++i) {
        Face& face = m_faces[i];
        
        if (face.edgeIndex >= edgeCount || face.edgeIndex == INVALID_INDEX) {
            // Face has an invalid edge reference
            // Try to locate a valid edge for this face
            for (size_t j = 0; j < edgeCount; ++j) {
                if (m_edges[j].leftFace == i || m_edges[j].rightFace == i) {
                    face.edgeIndex = static_cast<uint32_t>(j);
                    break;
                }
            }
            continue;  // Skip rest of checks if we couldn't find a valid edge
        }
        
        // Check face loop connectivity
        uint32_t firstEdge = face.edgeIndex;
        uint32_t currentEdge = firstEdge;
        bool isLeftFace = (m_edges[currentEdge].leftFace == i);
        
        std::vector<uint32_t> visitedEdges;
        std::vector<uint32_t> faceEdges;
        bool loopIsClosed = false;
        
        // Follow the edge loop and see if it returns to the first edge
        int safetyCounter = 0;
        do {
            visitedEdges.push_back(currentEdge);
            faceEdges.push_back(currentEdge);
            
            Edge& edge = m_edges[currentEdge];
            uint32_t nextEdge;
            
            if (isLeftFace) {
                nextEdge = edge.leftCCW;
            } else {
                nextEdge = edge.rightCW;
            }
            
            // Check for end of loop
            if (nextEdge == firstEdge) {
                loopIsClosed = true;
                break;
            }
            
            // Check for invalid next edge
            if (nextEdge == INVALID_INDEX || nextEdge >= edgeCount) {
                break;
            }
            
            // Check for cycle (not returning to first edge)
            if (std::find(visitedEdges.begin(), visitedEdges.end(), nextEdge) != visitedEdges.end()) {
                break;
            }
            
            currentEdge = nextEdge;
            safetyCounter++;
        } while (safetyCounter < 10); // Safety limit
        
        // If the loop is not closed, try to fix it
        if (!loopIsClosed && faceEdges.size() > 1) {
            // Connect the last edge to the first edge
            uint32_t lastEdge = faceEdges.back();
            
            if (isLeftFace) {
                m_edges[lastEdge].leftCCW = firstEdge;
                m_edges[firstEdge].leftCW = lastEdge;
            } else {
                m_edges[lastEdge].rightCW = firstEdge;
                m_edges[firstEdge].rightCCW = lastEdge;
            }
        }
    }
    
    // 3. Rebuild edge mapping (for faster edge finding)
    m_edgeMap.clear();
    for (size_t i = 0; i < edgeCount; ++i) {
        const Edge& edge = m_edges[i];
        m_edgeMap[edgeKey(edge.startVertex, edge.endVertex)] = static_cast<uint32_t>(i);
    }
    
    // 4. Verify each vertex has a valid edge reference
    for (size_t i = 0; i < vertexCount; ++i) {
        Vertex& vertex = m_vertices[i];
        
        if (vertex.edgeIndex >= edgeCount || vertex.edgeIndex == INVALID_INDEX) {
            // Vertex has an invalid edge reference
            // Find an edge that uses this vertex
            for (size_t j = 0; j < edgeCount; ++j) {
                if (m_edges[j].startVertex == i || m_edges[j].endVertex == i) {
                    vertex.edgeIndex = static_cast<uint32_t>(j);
                    break;
                }
            }
        }
    }
}

// Find edges connected to a specific vertex
std::vector<uint32_t> WingedEdgeMesh::findVertexEdges(uint32_t vertexIndex) const
{
    std::vector<uint32_t> connectedEdges;
    
    // Search through all edges for those that include this vertex
    for (uint32_t i = 0; i < m_edges.size(); ++i) {
        const Edge& edge = m_edges[i];
        if (edge.startVertex == vertexIndex || edge.endVertex == vertexIndex) {
            connectedEdges.push_back(i);
        }
    }
    
    return connectedEdges;
}
