#include "TriangleMesh.h"
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <cassert>
#include "geometry/vectors/intersections.h"
#include "Edges.h"

// Constructor
TriangleMesh::TriangleMesh() 
    : m_pVertexMesh(std::make_shared<Vertices>())
{
}

// Constructor with existing Vertices
TriangleMesh::TriangleMesh(std::shared_ptr<Vertices> vertexMesh)
    : m_pVertexMesh(vertexMesh)
{
    if (!m_pVertexMesh) {
        m_pVertexMesh = std::make_shared<Vertices>();
    }
}

// Bounding box
box3 TriangleMesh::getBox() const {
    return m_pVertexMesh->getBox();
}

// Clear all mesh data (vertices and triangles)
void TriangleMesh::clear() {
    m_pVertexMesh->clear();
    m_triangles.clear();
    invalidateEdges();
    incrementVersion();
}

// Get all vertices of a triangle
uint3 TriangleMesh::getTriangleVertices(uint32_t triangleIndex) const {
    return m_triangles[triangleIndex];
}

// Get number of triangles
uint32_t TriangleMesh::getTriangleCount() const {
    return static_cast<uint32_t>(m_triangles.size());
}

// Calculate the area of a triangle
double TriangleMesh::calculateTriangleArea(uint32_t triangleIndex) const {
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Convert float3 to double3 for precise calculations
    const double3 p1 = double3(m_pVertexMesh->getVertexPosition(verts.x));
    const double3 p2 = double3(m_pVertexMesh->getVertexPosition(verts.y));
    const double3 p3 = double3(m_pVertexMesh->getVertexPosition(verts.z));
    
    // Area = 0.5 * |cross(v1, v2)|
    return 0.5 * length(cross(p2 - p1, p3 - p1));
}

// Calculate triangle normal
double3 TriangleMesh::calculateTriangleNormal(uint32_t triangleIndex) const {
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Convert float3 to double3 for precise calculations
    const double3 p1 = double3(m_pVertexMesh->getVertexPosition(verts.x));
    const double3 p2 = double3(m_pVertexMesh->getVertexPosition(verts.y));
    const double3 p3 = double3(m_pVertexMesh->getVertexPosition(verts.z));
    
    // Normal = normalize(cross(v1, v2))
    double3 normal = cross(p2 - p1, p3 - p1);
    double len = length(normal);
    
    if (len > 1e-10) {
        return normal / len;
    } else {
        return double3(0.0, 0.0, 1.0); // Default normal if degenerate triangle
    }
}

// Compute barycentric coordinates for a point with respect to a triangle
float3 TriangleMesh::computeBary(uint32_t triangleIndex, const float3& point) const {
    const uint3 tri = getTriangleVertices(triangleIndex);
    const float3 v0 = m_pVertexMesh->getVertexPosition(tri.x);
    const float3 v1 = m_pVertexMesh->getVertexPosition(tri.y);
    const float3 v2 = m_pVertexMesh->getVertexPosition(tri.z);
    return computeBarycentricCoordinates(point, v0, v1, v2);
}

// Helper method for derived classes to add triangles directly to storage
uint32_t TriangleMesh::addTriangle(uint32_t v1, uint32_t v2, uint32_t v3) {
    m_triangles.emplace_back(v1, v2, v3);
    invalidateEdges();
    incrementVersion();
    return static_cast<uint32_t>(m_triangles.size() - 1);
}

// Extract triangles (move out, leaving vertices intact)
std::vector<uint3> TriangleMesh::extractTriangles() {
    std::vector<uint3> extracted = std::move(m_triangles);
    m_triangles.clear(); // Ensure triangles is in a valid empty state
    invalidateEdges();
    incrementVersion();
    return extracted;
}

// Static factory method to create an icosahedron
std::shared_ptr<TriangleMesh> TriangleMesh::createIcosahedron(double radius) {
    std::shared_ptr<TriangleMesh> mesh = std::make_shared<TriangleMesh>();
    
    // Golden ratio for icosahedron calculations
    const double PHI = 1.61803398874989484820;
    
    // Calculate vertex positions for a unit icosahedron
    double norm = std::sqrt(1.0 + PHI * PHI);
    float a = (float)(radius / norm);
    float b = (float)(radius * PHI / norm);
    
    // Add 12 vertices of the icosahedron
    mesh->m_pVertexMesh->addVertex(float3(0, a, b));  // 0
    mesh->m_pVertexMesh->addVertex(float3(0, a, -b)); // 1
    mesh->m_pVertexMesh->addVertex(float3(0, -a, b)); // 2
    mesh->m_pVertexMesh->addVertex(float3(0, -a, -b)); // 3
    mesh->m_pVertexMesh->addVertex(float3(a, b, 0));  // 4
    mesh->m_pVertexMesh->addVertex(float3(-a, b, 0)); // 5
    mesh->m_pVertexMesh->addVertex(float3(a, -b, 0)); // 6
    mesh->m_pVertexMesh->addVertex(float3(-a, -b, 0)); // 7
    mesh->m_pVertexMesh->addVertex(float3(b, 0, a));  // 8
    mesh->m_pVertexMesh->addVertex(float3(-b, 0, a)); // 9
    mesh->m_pVertexMesh->addVertex(float3(b, 0, -a)); // 10
    mesh->m_pVertexMesh->addVertex(float3(-b, 0, -a)); // 11
    
    // Add 20 triangular faces
    mesh->addTriangle(0, 8, 4);
    mesh->addTriangle(0, 4, 5);
    mesh->addTriangle(0, 5, 9);
    mesh->addTriangle(0, 9, 2);
    mesh->addTriangle(0, 2, 8);
    
    mesh->addTriangle(1, 5, 4);
    mesh->addTriangle(1, 4, 10);
    mesh->addTriangle(1, 10, 3);
    mesh->addTriangle(1, 3, 11);
    mesh->addTriangle(1, 11, 5);
    
    mesh->addTriangle(2, 7, 6);
    mesh->addTriangle(2, 6, 8);
    mesh->addTriangle(2, 9, 7);
    
    mesh->addTriangle(3, 6, 7);
    mesh->addTriangle(3, 7, 11);
    mesh->addTriangle(3, 10, 6);
    
    mesh->addTriangle(4, 8, 10);
    mesh->addTriangle(5, 11, 9);
    mesh->addTriangle(6, 10, 8);
    mesh->addTriangle(7, 9, 11);
    
    return mesh;
}

// Static factory method to create a sphere mesh
std::shared_ptr<TriangleMesh> TriangleMesh::createSphere(double radius, uint32_t subdivisionLevel) {
    // Create initial icosahedron as a TriangleMesh
    std::shared_ptr<TriangleMesh> triangleMesh = TriangleMesh::createIcosahedron(radius);
    
    triangleMesh->verifyTopology();
    
    // Subdivide the mesh the requested number of times
    for (uint32_t level = 0; level < subdivisionLevel; ++level) {
        triangleMesh = triangleMesh->subdivide();
        triangleMesh->verifyTopology();
    }
    
    triangleMesh->getOrCreateEdges();

    // Verify optimal tessellation
    triangleMesh->verifyTopology();
    
    // Additional verification for subdivided icosahedron: expected face count
    // An icosahedron has 20 faces; each subdivision multiplies faces by 4
    uint32_t F = triangleMesh->getTriangleCount();
    uint32_t expectedF = 20 * (1u << (2 * subdivisionLevel)); // 20 * 4^subdivisionLevel
    assert(F == expectedF && "Triangle count mismatch for subdivided icosahedron");
    
    return triangleMesh;
}

// Lazily compute and return edges
std::shared_ptr<Edges> TriangleMesh::getOrCreateEdges() {
    if (!m_pEdges) {
        m_pEdges = Edges::computeEdges(*this);
    }
    return m_pEdges;
}

std::shared_ptr<const Edges> TriangleMesh::getOrCreateEdges() const {
    if (!m_pEdges) {
        m_pEdges = Edges::computeEdges(*this);
    }
    return m_pEdges;
}

// Verify mesh topology using Euler's formula
void TriangleMesh::verifyTopology() const {
    // For a closed triangle mesh, Euler's formula: V - E + F = 2
    // We can verify both with and without edge information
    uint32_t V = m_pVertexMesh->getVertexCount();
    uint32_t F = getTriangleCount();
    
    // Basic verification: V = 2 + F/2 (derived from Euler's formula with E = 3F/2)
    assert(F % 2 == 0 && "Face count must be even for closed triangle mesh");
    uint32_t expectedV = 2 + F / 2;
    assert(V == expectedV && "Vertex count should equal 2 + F/2 for closed triangle mesh");
    
    // If edges are computed, verify them too
    if (m_pEdges) {
        uint32_t E = m_pEdges->getEdgeCount();
        uint32_t expectedE = (3 * F) / 2;
        assert(E == expectedE && "Edge count should equal 3F/2 for closed triangle mesh");
        assert(V - E + F == 2 && "Euler's formula V - E + F = 2 violated");
    }
}

// Subdivide the mesh once (creates a new subdivided mesh)
std::shared_ptr<TriangleMesh> TriangleMesh::subdivide() const {
    auto subdivided = std::make_shared<TriangleMesh>(m_pVertexMesh);
    
    // Map to store midpoints to avoid duplication
    std::unordered_map<uint64_t, uint32_t> midpoints;
    
    // Calculate average radius for vertex projection onto sphere
    double radius = 0.0;
    for (uint32_t i = 0; i < m_pVertexMesh->getVertexCount(); ++i) {
        radius += length(m_pVertexMesh->getVertexPosition(i));
    }
    radius /= m_pVertexMesh->getVertexCount();
    float fRadius = (float)radius;
    
    // Process each original triangle
    for (uint32_t triangleIdx = 0; triangleIdx < m_triangles.size(); ++triangleIdx) {
        const uint3& triangleVerts = m_triangles[triangleIdx];
        
        uint32_t v1 = triangleVerts.x;
        uint32_t v2 = triangleVerts.y;
        uint32_t v3 = triangleVerts.z;
        
        // Create midpoints
        uint32_t m12 = subdivided->getMidpoint(v1, v2, midpoints, fRadius);
        uint32_t m23 = subdivided->getMidpoint(v2, v3, midpoints, fRadius);
        uint32_t m31 = subdivided->getMidpoint(v3, v1, midpoints, fRadius);
        
        // Create four new triangular faces
        subdivided->addTriangle(v1, m12, m31);
        subdivided->addTriangle(m12, v2, m23);
        subdivided->addTriangle(m31, m23, v3);
        subdivided->addTriangle(m12, m23, m31);
    }
    
    return subdivided;
}

// Helper to get or create midpoint between two vertices
uint32_t TriangleMesh::getMidpoint(uint32_t v1, uint32_t v2, 
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
    const float3 pos1 = m_pVertexMesh->getVertexPosition(v1);
    const float3 pos2 = m_pVertexMesh->getVertexPosition(v2);
    
    // Calculate midpoint and project onto sphere
    float3 midpoint = (pos1 + pos2) * 0.5f;
    float len = length(midpoint);
    if (len > 1e-10) {
        midpoint = (midpoint / len) * radius;
    }
    
    // Add new vertex and register in midpoint map
    uint32_t midpointIdx = m_pVertexMesh->addVertex(midpoint);
    midpoints[key] = midpointIdx;
    
    return midpointIdx;
}

