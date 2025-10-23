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
    // Factory methods
    static std::shared_ptr<EdgeMesh> createSphere(double radius, uint32_t subdivisionLevel);

    // inherited TriangleMesh methods
    virtual void clear() override;
    virtual uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3) override;
    virtual std::vector<uint3> extractTriangles() override;

    // Edge access methods
    uint32_t getEdgeCount() const;
    std::pair<uint32_t, uint32_t> getEdge(uint32_t edgeIndex) const;
    
    // Topology verification
    virtual void verifyTopology() const override;

private:
    EdgeMesh();
    
    std::shared_ptr<Edges> m_pEdges;
};
