#include "TriangleMesh.h"
#include <algorithm>
#include <cmath>
#include "geometry/vectors/intersections.h"

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
    incrementVersion();
    return static_cast<uint32_t>(m_triangles.size() - 1);
}

// Extract triangles (move out, leaving vertices intact)
std::vector<uint3> TriangleMesh::extractTriangles() {
    std::vector<uint3> extracted = std::move(m_triangles);
    m_triangles.clear(); // Ensure triangles is in a valid empty state
    incrementVersion();
    return extracted;
}

// Static factory method to create an icosahedron
std::shared_ptr<TriangleMesh> TriangleMesh::createIcosahedron(double radius) {
    auto mesh = std::make_shared<TriangleMesh>();
    mesh->populateIcosahedron(radius);
    return mesh;
}

// Populate the mesh with icosahedron geometry
void TriangleMesh::populateIcosahedron(double radius) {
    clear();
    
    // Golden ratio for icosahedron calculations
    const double PHI = 1.61803398874989484820;
    
    // Calculate vertex positions for a unit icosahedron
    double norm = std::sqrt(1.0 + PHI * PHI);
    float a = (float)(radius / norm);
    float b = (float)(radius * PHI / norm);
    
    // Add 12 vertices of the icosahedron
    m_pVertexMesh->addVertex(float3(0, a, b));  // 0
    m_pVertexMesh->addVertex(float3(0, a, -b)); // 1
    m_pVertexMesh->addVertex(float3(0, -a, b)); // 2
    m_pVertexMesh->addVertex(float3(0, -a, -b)); // 3
    m_pVertexMesh->addVertex(float3(a, b, 0));  // 4
    m_pVertexMesh->addVertex(float3(-a, b, 0)); // 5
    m_pVertexMesh->addVertex(float3(a, -b, 0)); // 6
    m_pVertexMesh->addVertex(float3(-a, -b, 0)); // 7
    m_pVertexMesh->addVertex(float3(b, 0, a));  // 8
    m_pVertexMesh->addVertex(float3(-b, 0, a)); // 9
    m_pVertexMesh->addVertex(float3(b, 0, -a)); // 10
    m_pVertexMesh->addVertex(float3(-b, 0, -a)); // 11
    
    // Add 20 triangular faces
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

