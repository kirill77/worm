#pragma once

#include <vector>
#include <cstdint>
#include <memory>
#include "geometry/vectors/vector.h"
#include "VertexMesh.h"

class Mesh : public VertexMesh
{
public:
    // Constructors and main methods
    Mesh();
    virtual ~Mesh() = default;
    
    // Triangle operations (basic access, no connectivity)
    uint3 getTriangleVertices(uint32_t triangleIndex) const;
    uint32_t getTriangleCount() const;
    double calculateTriangleArea(uint32_t triangleIndex) const;
    double3 calculateTriangleNormal(uint32_t triangleIndex) const;
    
    // Compute barycentric coordinates of a point with respect to a triangle
    // Returns barycentric coordinates (w0, w1, w2) where point = w0*v0 + w1*v1 + w2*v2
    float3 computeBary(uint32_t triangleIndex, const float3& point) const;
    
    virtual uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3);

    // Clear mesh data (vertices and triangles)
    virtual void clear() override;
    
    // Extract triangles (move out, leaving vertices intact)
    virtual std::vector<uint3> extractTriangles();

protected:
    std::vector<uint3> m_triangles;
}; 