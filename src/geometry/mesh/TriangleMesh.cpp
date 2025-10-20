#include "TriangleMesh.h"
#include <algorithm>
#include <cmath>
#include "geometry/vectors/intersections.h"

// Constructor
TriangleMesh::TriangleMesh() {
}

// Clear all mesh data (vertices and triangles)
void TriangleMesh::clear() {
    VertexMesh::clear(); // Clear vertices (already increments version)
    m_triangles.clear();
    // Note: version already incremented by VertexMesh::clear()
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
    const double3 p1 = double3(getVertexPosition(verts.x));
    const double3 p2 = double3(getVertexPosition(verts.y));
    const double3 p3 = double3(getVertexPosition(verts.z));
    
    // Area = 0.5 * |cross(v1, v2)|
    return 0.5 * length(cross(p2 - p1, p3 - p1));
}

// Calculate triangle normal
double3 TriangleMesh::calculateTriangleNormal(uint32_t triangleIndex) const {
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Convert float3 to double3 for precise calculations
    const double3 p1 = double3(getVertexPosition(verts.x));
    const double3 p2 = double3(getVertexPosition(verts.y));
    const double3 p3 = double3(getVertexPosition(verts.z));
    
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
    const float3 v0 = getVertexPosition(tri.x);
    const float3 v1 = getVertexPosition(tri.y);
    const float3 v2 = getVertexPosition(tri.z);
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

