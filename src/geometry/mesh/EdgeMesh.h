#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <memory>
#include "geometry/vectors/vector.h"
#include "TriangleMesh.h"
#include "Edges.h"

class EdgeMesh : public TriangleMesh {
public:
    // Constructors and main methods
    EdgeMesh();
    EdgeMesh(double radius, uint32_t subdivisionLevel);

    // inherited TriangleMesh methods
    virtual void clear() override;
    virtual uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3) override;
    virtual std::vector<uint3> extractTriangles() override;

    // Edge access methods
    uint32_t getEdgeCount() const;
    std::pair<uint32_t, uint32_t> getEdge(uint32_t edgeIndex) const;

private:
    std::shared_ptr<Edges> m_pEdges;
    
    // Helper methods for icosahedron creation
    void createIcosahedron(double radius);
    void subdivide(uint32_t levels);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2, 
                       std::unordered_map<uint64_t, uint32_t>& midpoints,
                       float radius);
};
