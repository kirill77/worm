#pragma once

#include <vector>
#include <cstdint>
#include <memory>
#include "geometry/vectors/vector.h"
#include "Identifiable.h"
#include "geometry/vectors/box.h"

class Mesh : public Identifiable
{
public:
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    struct Vertex {
        float3 position;
        Vertex(const float3& pos) : position(pos) {}
    };

    // Constructors and main methods
    Mesh();
    virtual ~Mesh() = default;
    
    // Vertex operations
    uint32_t addVertex(const float3& position);
    float3 getVertexPosition(uint32_t index) const;
    void setVertexPosition(uint32_t index, const float3& position);
    uint32_t getVertexCount() const;
    
    // Triangle operations (basic access, no connectivity)
    uint3 getTriangleVertices(uint32_t triangleIndex) const;
    uint32_t getTriangleCount() const;
    double calculateTriangleArea(uint32_t triangleIndex) const;
    double3 calculateTriangleNormal(uint32_t triangleIndex) const;
    
    virtual uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3);

    // Clear mesh data
    virtual void clear();
    
    // Extract triangles (move out, leaving vertices intact)
    virtual std::vector<uint3> extractTriangles();

    // Version tracking
    uint64_t getVersion() const { return m_version; }

    // Bounding box (cached based on version)
    box3 getBox() const;

private:
    std::vector<Vertex> vertices;
    std::vector<uint3> triangles;
    uint64_t m_version = 0;
    
    // Cached bounding box
    mutable box3 m_cachedBox = box3::empty();
    mutable uint64_t m_cachedBoxVersion = UINT64_MAX; // Invalid version initially
}; 