#pragma once

#include <vector>
#include <unordered_map>
#include <memory>
#include <cstdint>

class TriangleMesh;

struct Edge {
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    uint32_t startVertex;
    uint32_t endVertex;
    
    Edge(uint32_t start, uint32_t end)
        : startVertex(start)
        , endVertex(end) {}
};

class Edges : public std::enable_shared_from_this<Edges> {
public:
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    static std::shared_ptr<Edges> computeEdges(const TriangleMesh& mesh);
    
    uint32_t getEdgeCount() const;
    std::pair<uint32_t, uint32_t> getEdge(uint32_t edgeIndex) const;
    
    static uint64_t directionlessEdgeKey(uint32_t v1, uint32_t v2);
    
private:
    Edges();
    
    std::vector<Edge> edges;
    
    uint32_t addEdge(uint32_t startVertex, uint32_t endVertex, std::unordered_map<uint64_t, uint32_t>& edgeMap);
    static uint64_t directionalEdgeKey(uint32_t startVertex, uint32_t endVertex);
};

