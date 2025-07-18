#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include "geometry/vectors/vector.h"
#include "Mesh.h"

class EdgeMesh : public Mesh {
public:
    struct Edge {
        uint32_t startVertex;
        uint32_t endVertex;
        uint32_t rightTriangle;
        uint32_t nextEdge;
        Edge(uint32_t start, uint32_t end) 
            : startVertex(start)
            , endVertex(end)
            , rightTriangle(INVALID_INDEX)
            , nextEdge(INVALID_INDEX) {}
    };

    // Constructors and main methods
    EdgeMesh();
    EdgeMesh(double radius, uint32_t subdivisionLevel);

    // inherited Mesh methods
    virtual void clear() override;
    virtual uint32_t addTriangle(uint32_t v1, uint32_t v2, uint32_t v3) override;
    virtual std::vector<uint3> extractTriangles() override;

    uint32_t addEdge(uint32_t startVertex, uint32_t endVertex);
    std::vector<uint32_t> getTriangleNeighbors(uint32_t triangleIndex) const;
    
    // Edge access methods
    uint32_t getEdgeCount() const;
    std::pair<uint32_t, uint32_t> getEdge(uint32_t edgeIndex) const;
    std::vector<std::pair<uint32_t, uint32_t>> getAllEdges() const;

private:
    std::vector<Edge> edges;
    std::unordered_map<uint64_t, uint32_t> edgeMap;

    static uint64_t directionalEdgeKey(uint32_t startVertex, uint32_t endVertex);
    static uint64_t directionlessEdgeKey(uint32_t v1, uint32_t v2);
    uint32_t findEdge(uint32_t startVertex, uint32_t endVertex) const;
    
    // Helper methods for icosahedron creation
    void createIcosahedron(double radius);
    void subdivide(uint32_t levels);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2, 
                       std::unordered_map<uint64_t, uint32_t>& midpoints,
                       double radius);
};
