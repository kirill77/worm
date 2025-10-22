#include "Edges.h"
#include "TriangleMesh.h"

Edges::Edges() {
}

std::shared_ptr<Edges> Edges::computeEdges(const TriangleMesh& mesh) {
    auto edges = std::shared_ptr<Edges>(new Edges());
    
    std::unordered_map<uint64_t, uint32_t> edgeMap;
    
    const uint32_t triangleCount = mesh.getTriangleCount();
    for (uint32_t i = 0; i < triangleCount; ++i) {
        uint3 triangle = mesh.getTriangleVertices(i);
        edges->addEdge(triangle.x, triangle.y, edgeMap);
        edges->addEdge(triangle.y, triangle.z, edgeMap);
        edges->addEdge(triangle.z, triangle.x, edgeMap);
    }
    
    return edges;
}

uint32_t Edges::addEdge(uint32_t startVertex, uint32_t endVertex, std::unordered_map<uint64_t, uint32_t>& edgeMap) {
    auto key = directionalEdgeKey(startVertex, endVertex);
    auto it = edgeMap.find(key);
    
    if (it != edgeMap.end()) {
        return it->second;
    }
    
    edges.emplace_back(startVertex, endVertex);
    uint32_t edgeIndex = static_cast<uint32_t>(edges.size() - 1);
    
    edgeMap[key] = edgeIndex;
    
    return edgeIndex;
}

uint32_t Edges::getEdgeCount() const {
    return static_cast<uint32_t>(edges.size());
}

std::pair<uint32_t, uint32_t> Edges::getEdge(uint32_t edgeIndex) const {
    if (edgeIndex < edges.size()) {
        const Edge& edge = edges[edgeIndex];
        return std::make_pair(edge.startVertex, edge.endVertex);
    }
    return std::make_pair(INVALID_INDEX, INVALID_INDEX);
}

uint64_t Edges::directionalEdgeKey(uint32_t startVertex, uint32_t endVertex) {
    return ((uint64_t)endVertex << 32) | (uint64_t)startVertex;
}

uint64_t Edges::directionlessEdgeKey(uint32_t v1, uint32_t v2) {
    return (v1 <= v2) ? directionalEdgeKey(v1, v2) : directionalEdgeKey(v2, v1);
}

