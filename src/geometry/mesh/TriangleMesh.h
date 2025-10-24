#pragma once

#include <vector>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include "geometry/vectors/vector.h"
#include "geometry/vectors/box.h"
#include "Vertices.h"
#include "Identifiable.h"

class Edges;

class TriangleMesh : public Identifiable
{
public:
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    // Constructors and main methods
    TriangleMesh();
    TriangleMesh(std::shared_ptr<Vertices> vertexMesh);
    virtual ~TriangleMesh() = default;
    
    // Factory methods
    static std::shared_ptr<TriangleMesh> createIcosahedron(double radius);
    static std::shared_ptr<TriangleMesh> createSphere(double radius, uint32_t subdivisionLevel);
    
    // Subdivision
    std::shared_ptr<TriangleMesh> subdivide() const;
    
    // Vertex mesh access
    std::shared_ptr<Vertices> getVertices() const { return m_pVertexMesh; }
    void setVertices(std::shared_ptr<Vertices> vertexMesh) { m_pVertexMesh = vertexMesh; }
    
    // Edge access (lazily computed)
    std::shared_ptr<Edges> getOrCreateEdges();
    std::shared_ptr<const Edges> getOrCreateEdges() const;
    
    // Triangle operations (basic access, no connectivity)
    uint3 getTriangleVertices(uint32_t triangleIndex) const;
    uint32_t getTriangleCount() const;
    double calculateTriangleArea(uint32_t triangleIndex) const;
    double3 calculateTriangleNormal(uint32_t triangleIndex) const;
    
    // Compute barycentric coordinates of a point with respect to a triangle
    // Returns barycentric coordinates (w0, w1, w2) where point = w0*v0 + w1*v1 + w2*v2
    float3 computeBary(uint32_t triangleIndex, const float3& point) const;
    
    uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3);

    // Clear mesh data (vertices and triangles)
    void clear();
    
    // Extract triangles (move out, leaving vertices intact)
    std::vector<uint3> extractTriangles();
    
    // Version tracking
    uint64_t getVersion() const { return m_version + m_pVertexMesh->getVersion(); }
    
    // Bounding box
    box3 getBox() const;
    
    // Topology verification
    void verifyTopology() const;

protected:
    void incrementVersion() { ++m_version; }
    void invalidateEdges() { m_pEdges.reset(); }
    uint32_t getMidpoint(uint32_t v1, uint32_t v2, 
                         std::unordered_map<uint64_t, uint32_t>& midpoints,
                         float radius);
    
    std::shared_ptr<Vertices> m_pVertexMesh;
    std::vector<uint3> m_triangles;
    mutable std::shared_ptr<Edges> m_pEdges;
    uint64_t m_version = 0;
};
