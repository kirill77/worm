#include "EdgeMesh.h"
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
    Mesh::clear(); // Clear base class data (vertices and triangles)
    edges.clear();
    edgeMap.clear();
}

// Get number of edges
uint32_t EdgeMesh::getEdgeCount() const {
    return static_cast<uint32_t>(edges.size());
}

// Get a specific edge as a pair of vertex indices
std::pair<uint32_t, uint32_t> EdgeMesh::getEdge(uint32_t edgeIndex) const {
    if (edgeIndex < edges.size()) {
        const Edge& edge = edges[edgeIndex];
        return std::make_pair(edge.startVertex, edge.endVertex);
    }
    return std::make_pair(INVALID_INDEX, INVALID_INDEX);
}

// Get all edges as pairs of vertex indices
std::vector<std::pair<uint32_t, uint32_t>> EdgeMesh::getAllEdges() const {
    std::vector<std::pair<uint32_t, uint32_t>> result;
    result.reserve(edges.size());
    
    for (const Edge& edge : edges) {
        result.emplace_back(edge.startVertex, edge.endVertex);
    }
    
    return result;
}

// Generate a key for edge lookup
uint64_t EdgeMesh::directionalEdgeKey(uint32_t startVertex, uint32_t endVertex) {
    return ((uint64_t)endVertex << 32) | (uint64_t)startVertex;
}
uint64_t EdgeMesh::directionlessEdgeKey(uint32_t v1, uint32_t v2) {
    return (v1 <= v2) ? directionalEdgeKey(v1, v2) : directionalEdgeKey(v2, v1);
}

// Find an edge by vertices
uint32_t EdgeMesh::findEdge(uint32_t startVertex, uint32_t endVertex) const {
    auto it = edgeMap.find(directionalEdgeKey(startVertex, endVertex));
    return (it != edgeMap.end()) ? it->second : INVALID_INDEX;
}

// Add an edge to the mesh
uint32_t EdgeMesh::addEdge(uint32_t startVertex, uint32_t endVertex) {
    auto key = directionalEdgeKey(startVertex, endVertex);
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

// Add a triangle to the mesh
uint32_t EdgeMesh::addTriangle(uint32_t v1, uint32_t v2, uint32_t v3)
{
#ifndef NDEBUG
    {
        const auto& p1 = getVertexPosition(v1);
        const auto& p2 = getVertexPosition(v2);
        const auto& p3 = getVertexPosition(v3);
        
        // Calculate triangle normal
        float3 normal = cross(p2 - p1, p3 - p1);
        
        // For a sphere mesh, the normal should point outward from the center
        // So the dot product with the position vector should be positive
        double dotProduct = dot(normal, p1);
        
        // Assert that the normal triangles outward (for convex meshes like spheres)
        assert(dotProduct > 0.0 && "Vertices must be in counter-clockwise order!");
    }
#endif
    
    // Store triangle vertices using base class method
    uint32_t triangleIndex = Mesh::addTriangle(v1, v2, v3);
    
    // Still maintain edges for neighbor queries and other functionality
    uint32_t e1 = addEdge(v1, v2);
    uint32_t e2 = addEdge(v2, v3);
    uint32_t e3 = addEdge(v3, v1);
    
    // Link edges to this triangle and to each other for neighbor traversal
    edges[e1].rightTriangle = triangleIndex;
    edges[e2].rightTriangle = triangleIndex;
    edges[e3].rightTriangle = triangleIndex;
    
    // Link edges in a loop: e1 -> e2 -> e3 -> e1
    edges[e1].nextEdge = e2;
    edges[e2].nextEdge = e3;
    edges[e3].nextEdge = e1;
    
    return triangleIndex;
}

// Extract triangles and clear edge connectivity data
std::vector<uint3> EdgeMesh::extractTriangles() {
    // Extract triangles using base class method
    std::vector<uint3> extracted = Mesh::extractTriangles();
    
    // Clear edge connectivity data since triangles are gone
    edges.clear();
    edgeMap.clear();
    
    return extracted;
}

// Get neighboring triangles
std::vector<uint32_t> EdgeMesh::getTriangleNeighbors(uint32_t triangleIndex) const {
    std::vector<uint32_t> neighbors;
    
    // Get triangle vertices
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Check each edge's opposite
    for (int i = 0; i < 3; ++i) {
        uint32_t v1 = verts[i];
        uint32_t v2 = verts[(i + 1) % 3];
        
        // Find edge going the other way
        uint32_t oppositeEdge = findEdge(v2, v1);
        if (oppositeEdge != INVALID_INDEX && oppositeEdge < edges.size()) {
            uint32_t neighborTriangle = edges[oppositeEdge].rightTriangle;
            if (neighborTriangle != INVALID_INDEX && neighborTriangle != triangleIndex) {
                // Avoid duplicates
                if (std::find(neighbors.begin(), neighbors.end(), neighborTriangle) == neighbors.end()) {
                    neighbors.push_back(neighborTriangle);
                }
            }
        }
    }
    
    return neighbors;
}



// Create an icosahedron with the given radius
void EdgeMesh::createIcosahedron(double radius) {
    clear();
    
    // Golden ratio for icosahedron calculations
    const double PHI = 1.61803398874989484820;
    
    // Calculate vertex positions for a unit icosahedron
    double norm = std::sqrt(1.0 + PHI * PHI);
    float a = (float)(radius / norm);
    float b = (float)(radius * PHI / norm);
    
    // Add 12 vertices of the icosahedron
    addVertex(float3(0, a, b));  // 0
    addVertex(float3(0, a, -b)); // 1
    addVertex(float3(0, -a, b)); // 2
    addVertex(float3(0, -a, -b)); // 3
    addVertex(float3(a, b, 0));  // 4
    addVertex(float3(-a, b, 0)); // 5
    addVertex(float3(a, -b, 0)); // 6
    addVertex(float3(-a, -b, 0)); // 7
    addVertex(float3(b, 0, a));  // 8
    addVertex(float3(-b, 0, a)); // 9
    addVertex(float3(b, 0, -a)); // 10
    addVertex(float3(-b, 0, -a)); // 11
    
    // Add 20 triangular triangles
    addTriangle(0, 8, 4);
    addTriangle(0, 4, 5);
    addTriangle(0, 5, 9);
    addTriangle(0, 9, 2);
    addTriangle(0, 2, 8);
    
    addTriangle(1, 5, 4);
    addTriangle(1, 4, 10);
    addTriangle(1, 10, 3);
    addTriangle(1, 3, 11);
    addTriangle(1, 11, 5);
    
    addTriangle(2, 7, 6);
    addTriangle(2, 6, 8);
    addTriangle(2, 9, 7);
    
    addTriangle(3, 6, 7);
    addTriangle(3, 7, 11);
    addTriangle(3, 10, 6);
    
    addTriangle(4, 8, 10);
    addTriangle(5, 11, 9);
    addTriangle(6, 10, 8);
    addTriangle(7, 9, 11);
}

// Subdivide the mesh
void EdgeMesh::subdivide(uint32_t levels) {
    if (levels == 0) return;
    
    for (uint32_t level = 0; level < levels; ++level) {
        // Extract original triangles (vertices remain in place)
        std::vector<uint3> originalTriangles = extractTriangles();
        
        // Map to store midpoints to avoid duplication
        std::unordered_map<uint64_t, uint32_t> midpoints;
        
        // Calculate average radius for vertex projection
        double radius = 0.0;
        for (uint32_t i = 0; i < getVertexCount(); ++i) {
            radius += length(getVertexPosition(i));
        }
        radius /= getVertexCount();
        
        // Clear edge data (triangles are already extracted)
        edges.clear();
        edgeMap.clear();
        
        // Process each original triangle
        for (uint32_t triangleIdx = 0; triangleIdx < originalTriangles.size(); ++triangleIdx) {
            const uint3& triangleVerts = originalTriangles[triangleIdx];
            
            uint32_t v1 = triangleVerts.x;
            uint32_t v2 = triangleVerts.y;
            uint32_t v3 = triangleVerts.z;
            
            // Create midpoints
            float fRadius = (float)radius;
            uint32_t m12 = getMidpoint(v1, v2, midpoints, fRadius);
            uint32_t m23 = getMidpoint(v2, v3, midpoints, fRadius);
            uint32_t m31 = getMidpoint(v3, v1, midpoints, fRadius);
            
            // Create four new triangular triangles
            addTriangle(v1, m12, m31);
            addTriangle(m12, v2, m23);
            addTriangle(m31, m23, v3);
            addTriangle(m12, m23, m31);
        }
    }
}

// Helper to get or create midpoint between two vertices
uint32_t EdgeMesh::getMidpoint(uint32_t v1, uint32_t v2, 
                               std::unordered_map<uint64_t, uint32_t>& midpoints,
                               float radius) {
    // Use a normalized edge key (smaller vertex index first)
    uint64_t key = directionlessEdgeKey(v1, v2);
    
    // Check if midpoint already exists
    auto it = midpoints.find(key);
    if (it != midpoints.end()) {
        return it->second;
    }
    
    // Create a new midpoint using current mesh vertex positions
    const auto &pos1 = getVertexPosition(v1);
    const auto &pos2 = getVertexPosition(v2);
    
    // Calculate midpoint and project onto sphere
    float3 midpoint = (pos1 + pos2) * 0.5f;
    float len = length(midpoint);
    if (len > 1e-10) {
        midpoint = (midpoint / len) * radius;
    }
    
    // Add new vertex and register in midpoint map
    uint32_t midpointIdx = addVertex(midpoint);
    midpoints[key] = midpointIdx;
    
    return midpointIdx;
}



