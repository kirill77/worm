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

// Static factory method to create a sphere mesh
std::shared_ptr<EdgeMesh> EdgeMesh::createSphere(double radius, uint32_t subdivisionLevel) {
    // Create initial icosahedron as a TriangleMesh
    std::shared_ptr<TriangleMesh> triangleMesh = TriangleMesh::createIcosahedron(radius);

    triangleMesh->verifyTopology();
    
    // Subdivide the mesh the requested number of times
    for (uint32_t level = 0; level < subdivisionLevel; ++level) {
        triangleMesh = triangleMesh->subdivide();
        triangleMesh->verifyTopology();
    }
    
    // Create EdgeMesh and transfer geometry from TriangleMesh
    std::shared_ptr<EdgeMesh> edgeMesh = std::shared_ptr<EdgeMesh>(new EdgeMesh());
    edgeMesh->m_triangles = triangleMesh->extractTriangles();
    edgeMesh->setVertices(triangleMesh->getVertices());
    edgeMesh->incrementVersion();
    
    // Compute edges from all triangles
    edgeMesh->m_pEdges = Edges::computeEdges(*edgeMesh);
    
    // Verify optimal tessellation
    edgeMesh->verifyTopology();
    
    // Additional verification for subdivided icosahedron: expected face count
    // An icosahedron has 20 faces; each subdivision multiplies faces by 4
    uint32_t F = edgeMesh->getTriangleCount();
    uint32_t expectedF = 20 * (1u << (2 * subdivisionLevel)); // 20 * 4^subdivisionLevel
    assert(F == expectedF && "Triangle count mismatch for subdivided icosahedron");
    
    return edgeMesh;
}

// Verify mesh topology using Euler's formula with edge information
void EdgeMesh::verifyTopology() const {
    // Call base class verification first
    TriangleMesh::verifyTopology();
    
    // Additional verification with edge count
    // For a closed triangle mesh: E = 3F/2 (each triangle has 3 edges, each edge shared by 2)
    // Euler's formula: V - E + F = 2
    uint32_t V = getVertices()->getVertexCount();
    uint32_t E = getEdgeCount();
    uint32_t F = getTriangleCount();
    uint32_t expectedE = (3 * F) / 2;
    
    assert(E == expectedE && "Edge count should equal 3F/2 for closed triangle mesh");
    assert(V - E + F == 2 && "Euler's formula V - E + F = 2 violated");
}



// Clear all mesh data
void EdgeMesh::clear() {
    TriangleMesh::clear(); // Clear base class data (vertices and triangles)
    m_pEdges.reset();
}

// Get number of edges
uint32_t EdgeMesh::getEdgeCount() const {
    return m_pEdges ? m_pEdges->getEdgeCount() : 0;
}

// Get a specific edge as a pair of vertex indices
std::pair<uint32_t, uint32_t> EdgeMesh::getEdge(uint32_t edgeIndex) const {
    if (m_pEdges) {
        return m_pEdges->getEdge(edgeIndex);
    }
    return std::make_pair(INVALID_INDEX, INVALID_INDEX);
}

// Add a triangle to the mesh
uint32_t EdgeMesh::addTriangle(uint32_t v1, uint32_t v2, uint32_t v3)
{
#ifndef NDEBUG
    {
        const auto& p1 = getVertices()->getVertexPosition(v1);
        const auto& p2 = getVertices()->getVertexPosition(v2);
        const auto& p3 = getVertices()->getVertexPosition(v3);
        
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
    
    return triangleIndex;
}

// Extract triangles and clear edge connectivity data
std::vector<uint3> EdgeMesh::extractTriangles() {
    // Extract triangles using base class method
    std::vector<uint3> extracted = TriangleMesh::extractTriangles();
    
    // Clear edge connectivity data since triangles are gone
    m_pEdges.reset();
    
    return extracted;
}
