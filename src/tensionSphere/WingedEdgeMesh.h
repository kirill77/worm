#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include "../math/vector.h"

class WingedEdgeMesh {
public:
    static const uint32_t INVALID_INDEX = UINT32_MAX;
    
    struct Vertex {
        double3 position;
        Vertex(const double3& pos) : position(pos) {}
    };

    struct Edge {
        uint32_t startVertex;
        uint32_t endVertex;
        uint32_t rightFace;
        uint32_t nextEdge;
        Edge(uint32_t start, uint32_t end) 
            : startVertex(start)
            , endVertex(end)
            , rightFace(INVALID_INDEX)
            , nextEdge(INVALID_INDEX) {}
    };

    struct Face {
        uint32_t edgeIndex;
        Face(uint32_t edge) : edgeIndex(edge) {}
    };

    // Constructors and main methods
    WingedEdgeMesh();
    WingedEdgeMesh(double radius, uint32_t subdivisionLevel);
    void clear();
    uint32_t addVertex(const double3& position);
    uint32_t addEdge(uint32_t startVertex, uint32_t endVertex);
    uint32_t addFace(uint32_t v1, uint32_t v2, uint32_t v3);
    std::vector<uint32_t> getFaceVertices(uint32_t faceIndex) const;
    std::vector<uint32_t> getFaceNeighbors(uint32_t faceIndex) const;
    double calculateFaceArea(uint32_t faceIndex) const;
    double3 calculateFaceNormal(uint32_t faceIndex) const;
    
    // Additional accessor methods
    double3 getVertexPosition(uint32_t index) const;
    void setVertexPosition(uint32_t index, const double3& position);
    uint32_t getVertexCount() const;
    uint32_t getFaceCount() const;

private:
    std::vector<Vertex> vertices;
    std::vector<Edge> edges;
    std::vector<Face> faces;
    std::unordered_map<std::string, uint32_t> edgeMap;

    std::string edgeKey(uint32_t startVertex, uint32_t endVertex) const;
    uint32_t findEdge(uint32_t startVertex, uint32_t endVertex) const;
    
    // Helper methods for icosahedron creation
    void createIcosahedron(double radius);
    void subdivide(uint32_t levels);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2, 
                       std::unordered_map<std::string, uint32_t>& midpoints,
                       double radius);
};
