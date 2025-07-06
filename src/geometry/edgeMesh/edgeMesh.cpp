#include "edgeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cassert>  // For assert()

// Constructor
EdgeMesh::EdgeMesh() {
}

// Constructor with radius and subdivision
EdgeMesh::EdgeMesh(double radius, uint32_t subdivisionLevel) {
    createIcosahedron(radius);
    if (subdivisionLevel > 0) {
        subdivide(subdivisionLevel);
    }
}

// Clear all mesh data
void EdgeMesh::clear() {
    vertices.clear();
    edges.clear();
    faces.clear();
    edgeMap.clear();
}

// Add a vertex to the mesh
uint32_t EdgeMesh::addVertex(const double3& position) {
    vertices.emplace_back(position);
    return static_cast<uint32_t>(vertices.size() - 1);
}

// Get vertex position
double3 EdgeMesh::getVertexPosition(uint32_t index) const {
    if (index < vertices.size()) {
        return vertices[index].position;
    }
    return double3(0.0, 0.0, 0.0); // Return zero vector for invalid index
}

// Set vertex position
void EdgeMesh::setVertexPosition(uint32_t index, const double3& position) {
    if (index < vertices.size()) {
        vertices[index].position = position;
    }
}

// Get number of vertices
uint32_t EdgeMesh::getVertexCount() const {
    return static_cast<uint32_t>(vertices.size());
}

// Get number of faces
uint32_t EdgeMesh::getFaceCount() const {
    return static_cast<uint32_t>(faces.size());
}

// Generate a key for edge lookup
uint64_t EdgeMesh::edgeKey(uint32_t startVertex, uint32_t endVertex) const {
    return ((uint64_t)endVertex << 32) | (uint64_t)startVertex;
}

// Find an edge by vertices
uint32_t EdgeMesh::findEdge(uint32_t startVertex, uint32_t endVertex) const {
    auto it = edgeMap.find(edgeKey(startVertex, endVertex));
    return (it != edgeMap.end()) ? it->second : INVALID_INDEX;
}

// Add an edge to the mesh
uint32_t EdgeMesh::addEdge(uint32_t startVertex, uint32_t endVertex) {
    auto key = edgeKey(startVertex, endVertex);
    auto it = edgeMap.find(key);
    
    if (it != edgeMap.end()) {
        return it->second; // Edge already exists
    }
    
    // Create new edge
    edges.emplace_back(startVertex, endVertex);
    uint32_t edgeIndex = static_cast<uint32_t>(edges.size() - 1);
    
    // Add to map for fast lookup
    edgeMap[key] = edgeIndex;
    
    return edgeIndex;
}

// Add a face to the mesh
uint32_t EdgeMesh::addFace(uint32_t v1, uint32_t v2, uint32_t v3)
{
#ifndef NDEBUG
    {
        const double3& p1 = vertices[v1].position;
        const double3& p2 = vertices[v2].position;
        const double3& p3 = vertices[v3].position;
        
        // Calculate face normal
        double3 normal = cross(p2 - p1, p3 - p1);
        
        // For a sphere mesh, the normal should point outward from the center
        // So the dot product with the position vector should be positive
        double dotProduct = dot(normal, p1);
        
        // Assert that the normal faces outward (for convex meshes like spheres)
        assert(dotProduct > 0.0 && "Vertices must be in counter-clockwise order!");
    }
#endif
    
    // Add or find edges
    uint32_t e1 = addEdge(v1, v2);
    uint32_t e2 = addEdge(v2, v3);
    uint32_t e3 = addEdge(v3, v1);
    
    // Create new face
    faces.emplace_back(e1);
    uint32_t faceIndex = static_cast<uint32_t>(faces.size() - 1);
    
    // Link edges to this face and to each other
    edges[e1].rightFace = faceIndex;
    edges[e2].rightFace = faceIndex;
    edges[e3].rightFace = faceIndex;
    
    // Link edges in a loop: e1 -> e2 -> e3 -> e1
    edges[e1].nextEdge = e2;
    edges[e2].nextEdge = e3;
    edges[e3].nextEdge = e1;
    
    return faceIndex;
}

// Get all vertices of a face
std::vector<uint32_t> EdgeMesh::getFaceVertices(uint32_t faceIndex) const {
    std::vector<uint32_t> result;
    
    if (faceIndex >= faces.size()) {
        return result;
    }
    
    uint32_t startEdgeIndex = faces[faceIndex].edgeIndex;
    if (startEdgeIndex >= edges.size()) {
        return result;
    }
    
    // Follow the edge loop and collect vertices
    uint32_t currentEdge = startEdgeIndex;
    for (int i = 0; i < 3; ++i) { // assuming triangular faces
        if (currentEdge == INVALID_INDEX || currentEdge >= edges.size()) {
            break;
        }
        
        result.push_back(edges[currentEdge].startVertex);
        currentEdge = edges[currentEdge].nextEdge;
    }
    
    return result;
}

// Get neighboring faces
std::vector<uint32_t> EdgeMesh::getFaceNeighbors(uint32_t faceIndex) const {
    std::vector<uint32_t> neighbors;
    
    if (faceIndex >= faces.size()) {
        return neighbors;
    }
    
    // Get face vertices
    std::vector<uint32_t> verts = getFaceVertices(faceIndex);
    if (verts.size() < 3) {
        return neighbors;
    }
    
    // Check each edge's opposite
    for (size_t i = 0; i < verts.size(); ++i) {
        uint32_t v1 = verts[i];
        uint32_t v2 = verts[(i + 1) % verts.size()];
        
        // Find edge going the other way
        uint32_t oppositeEdge = findEdge(v2, v1);
        if (oppositeEdge != INVALID_INDEX && oppositeEdge < edges.size()) {
            uint32_t neighborFace = edges[oppositeEdge].rightFace;
            if (neighborFace != INVALID_INDEX && neighborFace != faceIndex) {
                // Avoid duplicates
                if (std::find(neighbors.begin(), neighbors.end(), neighborFace) == neighbors.end()) {
                    neighbors.push_back(neighborFace);
                }
            }
        }
    }
    
    return neighbors;
}

// Calculate the area of a face
double EdgeMesh::calculateFaceArea(uint32_t faceIndex) const {
    std::vector<uint32_t> verts = getFaceVertices(faceIndex);
    if (verts.size() < 3) {
        return 0.0;
    }
    
    const double3& p1 = vertices[verts[0]].position;
    const double3& p2 = vertices[verts[1]].position;
    const double3& p3 = vertices[verts[2]].position;
    
    // Area = 0.5 * |cross(v1, v2)|
    return 0.5 * length(cross(p2 - p1, p3 - p1));
}

// Calculate face normal
double3 EdgeMesh::calculateFaceNormal(uint32_t faceIndex) const {
    std::vector<uint32_t> verts = getFaceVertices(faceIndex);
    if (verts.size() < 3) {
        return double3(0.0, 0.0, 1.0); // Default normal if face is invalid
    }
    
    const double3& p1 = vertices[verts[0]].position;
    const double3& p2 = vertices[verts[1]].position;
    const double3& p3 = vertices[verts[2]].position;
    
    // Normal = normalize(cross(v1, v2))
    double3 normal = cross(p2 - p1, p3 - p1);
    double len = length(normal);
    
    if (len > 1e-10) {
        return normal / len;
    } else {
        return double3(0.0, 0.0, 1.0); // Default normal if degenerate face
    }
}

// Create an icosahedron with the given radius
void EdgeMesh::createIcosahedron(double radius) {
    clear();
    
    // Golden ratio for icosahedron calculations
    const double PHI = 1.61803398874989484820;
    
    // Calculate vertex positions for a unit icosahedron
    double norm = std::sqrt(1.0 + PHI * PHI);
    double a = radius / norm;
    double b = radius * PHI / norm;
    
    // Add 12 vertices of the icosahedron
    addVertex(double3(0, a, b));  // 0
    addVertex(double3(0, a, -b)); // 1
    addVertex(double3(0, -a, b)); // 2
    addVertex(double3(0, -a, -b)); // 3
    addVertex(double3(a, b, 0));  // 4
    addVertex(double3(-a, b, 0)); // 5
    addVertex(double3(a, -b, 0)); // 6
    addVertex(double3(-a, -b, 0)); // 7
    addVertex(double3(b, 0, a));  // 8
    addVertex(double3(-b, 0, a)); // 9
    addVertex(double3(b, 0, -a)); // 10
    addVertex(double3(-b, 0, -a)); // 11
    
    // Add 20 triangular faces
    addFace(0, 8, 4);
    addFace(0, 4, 5);
    addFace(0, 5, 9);
    addFace(0, 9, 2);
    addFace(0, 2, 8);
    
    addFace(1, 5, 4);
    addFace(1, 4, 10);
    addFace(1, 10, 3);
    addFace(1, 3, 11);
    addFace(1, 11, 5);
    
    addFace(2, 7, 6);
    addFace(2, 6, 8);
    addFace(2, 9, 7);
    
    addFace(3, 6, 7);
    addFace(3, 7, 11);
    addFace(3, 10, 6);
    
    addFace(4, 8, 10);
    addFace(5, 11, 9);
    addFace(6, 10, 8);
    addFace(7, 9, 11);
}

// Subdivide the mesh
void EdgeMesh::subdivide(uint32_t levels) {
    if (levels == 0) return;
    
    for (uint32_t level = 0; level < levels; ++level) {
        // Store original mesh data
        std::vector<Face> originalFaces = faces;
        std::vector<Vertex> originalVertices = vertices;
        
        // Map to store midpoints to avoid duplication
        std::unordered_map<std::string, uint32_t> midpoints;
        
        // Calculate average radius for vertex projection
        double radius = 0.0;
        for (const auto& vertex : vertices) {
            radius += length(vertex.position);
        }
        radius /= vertices.size();
        
        // Get face vertex indices before clearing
        std::vector<std::vector<uint32_t>> faceVertices;
        for (uint32_t i = 0; i < originalFaces.size(); ++i) {
            faceVertices.push_back(getFaceVertices(i));
        }
        
        // Clear faces and edge data
        faces.clear();
        edges.clear();
        edgeMap.clear();
        
        // Process each original face
        for (size_t faceIdx = 0; faceIdx < faceVertices.size(); ++faceIdx) {
            const auto& faceVerts = faceVertices[faceIdx];
            
            if (faceVerts.size() != 3) {
                // Skip invalid faces
                continue;
            }
            
            uint32_t v1 = faceVerts[0];
            uint32_t v2 = faceVerts[1];
            uint32_t v3 = faceVerts[2];
            
            // Create midpoints
            uint32_t m12 = getMidpoint(v1, v2, midpoints, radius);
            uint32_t m23 = getMidpoint(v2, v3, midpoints, radius);
            uint32_t m31 = getMidpoint(v3, v1, midpoints, radius);
            
            // Create four new triangular faces
            addFace(v1, m12, m31);
            addFace(m12, v2, m23);
            addFace(m31, m23, v3);
            addFace(m12, m23, m31);
        }
    }
}

// Helper to get or create midpoint between two vertices
uint32_t EdgeMesh::getMidpoint(uint32_t v1, uint32_t v2, 
                               std::unordered_map<std::string, uint32_t>& midpoints,
                               double radius) {
    // Use a normalized edge key (smaller vertex index first)
    std::string key = (v1 < v2) ? 
        (std::to_string(v1) + ":" + std::to_string(v2)) : 
        (std::to_string(v2) + ":" + std::to_string(v1));
    
    // Check if midpoint already exists
    auto it = midpoints.find(key);
    if (it != midpoints.end()) {
        return it->second;
    }
    
    // Create a new midpoint
    double3 pos1 = vertices[v1].position;
    double3 pos2 = vertices[v2].position;
    
    // Calculate midpoint and project onto sphere
    double3 midpoint = (pos1 + pos2) * 0.5;
    double len = length(midpoint);
    if (len > 1e-10) {
        midpoint = (midpoint / len) * radius;
    }
    
    // Add new vertex and register in midpoint map
    uint32_t midpointIdx = addVertex(midpoint);
    midpoints[key] = midpointIdx;
    
    return midpointIdx;
}

