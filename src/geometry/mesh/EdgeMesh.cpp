#include "EdgeMesh.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cassert>  // For assert()

// Constructor
EdgeMesh::EdgeMesh()
    : m_pEdges(std::make_shared<Edges>()) {
}

// Constructor with radius and subdivision
EdgeMesh::EdgeMesh(double radius, uint32_t subdivisionLevel)
    : m_pEdges(std::make_shared<Edges>()) {
    createIcosahedron(radius);
    if (subdivisionLevel > 0) {
        subdivide(subdivisionLevel);
    }
}



// Clear all mesh data
void EdgeMesh::clear() {
    TriangleMesh::clear(); // Clear base class data (vertices and triangles)
    m_pEdges->clear();
}

// Get number of edges
uint32_t EdgeMesh::getEdgeCount() const {
    return m_pEdges->getEdgeCount();
}

// Get a specific edge as a pair of vertex indices
std::pair<uint32_t, uint32_t> EdgeMesh::getEdge(uint32_t edgeIndex) const {
    return m_pEdges->getEdge(edgeIndex);
}

// Add an edge to the mesh
uint32_t EdgeMesh::addEdge(uint32_t startVertex, uint32_t endVertex) {
    return m_pEdges->addEdge(startVertex, endVertex);
}

// Add a triangle to the mesh
uint32_t EdgeMesh::addTriangle(uint32_t v1, uint32_t v2, uint32_t v3)
{
#ifndef NDEBUG
    {
        const auto& p1 = getVertexMesh()->getVertexPosition(v1);
        const auto& p2 = getVertexMesh()->getVertexPosition(v2);
        const auto& p3 = getVertexMesh()->getVertexPosition(v3);
        
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
    uint32_t triangleIndex = TriangleMesh::addTriangle(v1, v2, v3);
    
    // Still maintain edges for neighbor queries and other functionality
    addEdge(v1, v2);
    addEdge(v2, v3);
    addEdge(v3, v1);
    
    return triangleIndex;
}

// Extract triangles and clear edge connectivity data
std::vector<uint3> EdgeMesh::extractTriangles() {
    // Extract triangles using base class method
    std::vector<uint3> extracted = TriangleMesh::extractTriangles();
    
    // Clear edge connectivity data since triangles are gone
    m_pEdges->clear();
    
    return extracted;
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
    getVertexMesh()->addVertex(float3(0, a, b));  // 0
    getVertexMesh()->addVertex(float3(0, a, -b)); // 1
    getVertexMesh()->addVertex(float3(0, -a, b)); // 2
    getVertexMesh()->addVertex(float3(0, -a, -b)); // 3
    getVertexMesh()->addVertex(float3(a, b, 0));  // 4
    getVertexMesh()->addVertex(float3(-a, b, 0)); // 5
    getVertexMesh()->addVertex(float3(a, -b, 0)); // 6
    getVertexMesh()->addVertex(float3(-a, -b, 0)); // 7
    getVertexMesh()->addVertex(float3(b, 0, a));  // 8
    getVertexMesh()->addVertex(float3(-b, 0, a)); // 9
    getVertexMesh()->addVertex(float3(b, 0, -a)); // 10
    getVertexMesh()->addVertex(float3(-b, 0, -a)); // 11
    
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
        for (uint32_t i = 0; i < getVertexMesh()->getVertexCount(); ++i) {
            radius += length(getVertexMesh()->getVertexPosition(i));
        }
        radius /= getVertexMesh()->getVertexCount();
        
        // Clear edge data (triangles are already extracted)
        m_pEdges->clear();
        
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
    uint64_t key = Edges::directionlessEdgeKey(v1, v2);
    
    // Check if midpoint already exists
    auto it = midpoints.find(key);
    if (it != midpoints.end()) {
        return it->second;
    }
    
    // Create a new midpoint using current mesh vertex positions
    const auto &pos1 = getVertexMesh()->getVertexPosition(v1);
    const auto &pos2 = getVertexMesh()->getVertexPosition(v2);
    
    // Calculate midpoint and project onto sphere
    float3 midpoint = (pos1 + pos2) * 0.5f;
    float len = length(midpoint);
    if (len > 1e-10) {
        midpoint = (midpoint / len) * radius;
    }
    
    // Add new vertex and register in midpoint map
    uint32_t midpointIdx = getVertexMesh()->addVertex(midpoint);
    midpoints[key] = midpointIdx;
    
    return midpointIdx;
}



