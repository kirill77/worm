#include "Vertices.h"
#include <algorithm>

// Constructor
Vertices::Vertices() {
}

// Get bounding box (cached based on version)
box3 Vertices::getBox() const {
    // Check if cached box is valid for current version
    if (m_cachedBoxVersion == m_version) {
        return m_cachedBox;
    }
    
    // Recompute bounding box
    if (m_vertices.empty()) {
        m_cachedBox = box3::empty();
    } else {
        // Initialize with first vertex
        float3 minPt = m_vertices[0].position;
        float3 maxPt = m_vertices[0].position;
        
        // Find min/max for all vertices
        for (size_t i = 1; i < m_vertices.size(); ++i) {
            const float3& pos = m_vertices[i].position;
            minPt = min(minPt, pos);
            maxPt = max(maxPt, pos);
        }
        
        m_cachedBox = box3(minPt, maxPt);
    }
    
    // Update cached version
    m_cachedBoxVersion = m_version;
    
    return m_cachedBox;
}

// Clear all vertex data
void Vertices::clear() {
    m_vertices.clear();
    ++m_version;
}

// Add a vertex to the mesh
uint32_t Vertices::addVertex(const float3& position) {
    m_vertices.emplace_back(position);
    ++m_version;
    return static_cast<uint32_t>(m_vertices.size() - 1);
}

// Get vertex position
float3 Vertices::getVertexPosition(uint32_t index) const {
    if (index < m_vertices.size()) {
        return m_vertices[index].position;
    }
    return float3(0.0f, 0.0f, 0.0f); // Return zero vector for invalid index
}

// Set vertex position
void Vertices::setVertexPosition(uint32_t index, const float3& position) {
    if (index < m_vertices.size()) {
        m_vertices[index].position = position;
        ++m_version;
    }
}

// Get number of vertices
uint32_t Vertices::getVertexCount() const {
    return static_cast<uint32_t>(m_vertices.size());
}

