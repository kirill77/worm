#include "Mesh.h"
#include <algorithm>
#include <cmath>

// Constructor
Mesh::Mesh() {
}

// Get bounding box (cached based on version)
box3 Mesh::getBox() const {
    // Check if cached box is valid for current version
    if (m_cachedBoxVersion == m_version) {
        return m_cachedBox;
    }
    
    // Recompute bounding box
    if (vertices.empty()) {
        m_cachedBox = box3::empty();
    } else {
        // Initialize with first vertex
        float3 minPt = vertices[0].position;
        float3 maxPt = vertices[0].position;
        
        // Find min/max for all vertices
        for (size_t i = 1; i < vertices.size(); ++i) {
            const float3& pos = vertices[i].position;
            minPt = min(minPt, pos);
            maxPt = max(maxPt, pos);
        }
        
        m_cachedBox = box3(minPt, maxPt);
    }
    
    // Update cached version
    m_cachedBoxVersion = m_version;
    
    return m_cachedBox;
}

// Clear all mesh data
void Mesh::clear() {
    vertices.clear();
    triangles.clear();
    ++m_version;
}

// Add a vertex to the mesh
uint32_t Mesh::addVertex(const float3& position) {
    vertices.emplace_back(position);
    ++m_version;
    return static_cast<uint32_t>(vertices.size() - 1);
}

// Get vertex position
float3 Mesh::getVertexPosition(uint32_t index) const {
    if (index < vertices.size()) {
        return vertices[index].position;
    }
    return float3(0.0f, 0.0f, 0.0f); // Return zero vector for invalid index
}

// Set vertex position
void Mesh::setVertexPosition(uint32_t index, const float3& position) {
    if (index < vertices.size()) {
        vertices[index].position = position;
        ++m_version;
    }
}

// Get number of vertices
uint32_t Mesh::getVertexCount() const {
    return static_cast<uint32_t>(vertices.size());
}

// Get all vertices of a triangle
uint3 Mesh::getTriangleVertices(uint32_t triangleIndex) const {
    return triangles[triangleIndex];
}

// Get number of triangles
uint32_t Mesh::getTriangleCount() const {
    return static_cast<uint32_t>(triangles.size());
}

// Calculate the area of a triangle
double Mesh::calculateTriangleArea(uint32_t triangleIndex) const {
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Convert float3 to double3 for precise calculations
    const double3 p1 = double3(vertices[verts.x].position);
    const double3 p2 = double3(vertices[verts.y].position);
    const double3 p3 = double3(vertices[verts.z].position);
    
    // Area = 0.5 * |cross(v1, v2)|
    return 0.5 * length(cross(p2 - p1, p3 - p1));
}

// Calculate triangle normal
double3 Mesh::calculateTriangleNormal(uint32_t triangleIndex) const {
    uint3 verts = getTriangleVertices(triangleIndex);
    
    // Convert float3 to double3 for precise calculations
    const double3 p1 = double3(vertices[verts.x].position);
    const double3 p2 = double3(vertices[verts.y].position);
    const double3 p3 = double3(vertices[verts.z].position);
    
    // Normal = normalize(cross(v1, v2))
    double3 normal = cross(p2 - p1, p3 - p1);
    double len = length(normal);
    
    if (len > 1e-10) {
        return normal / len;
    } else {
        return double3(0.0, 0.0, 1.0); // Default normal if degenerate triangle
    }
}

// Helper method for derived classes to add triangles directly to storage
uint32_t Mesh::addTriangle(uint32_t v1, uint32_t v2, uint32_t v3) {
    triangles.emplace_back(v1, v2, v3);
    ++m_version;
    return static_cast<uint32_t>(triangles.size() - 1);
}

// Extract triangles (move out, leaving vertices intact)
std::vector<uint3> Mesh::extractTriangles() {
    std::vector<uint3> extracted = std::move(triangles);
    triangles.clear(); // Ensure triangles is in a valid empty state
    ++m_version;
    return extracted;
}



 