#pragma once

#include <vector>
#include <unordered_map>
#include <string>
#include <cstdint>
#include "../math/vector.h"

// A class that implements a Winged-Edge data structure for mesh representation
class WingedEdgeMesh {
public:
    // Forward declarations for Winged-Edge data structure components
    struct Vertex;
    struct Edge;
    struct Face;
    
    // Constructor creates an empty mesh
    WingedEdgeMesh();
    
    // Constructor creates a mesh from an icosahedron with optional subdivision
    WingedEdgeMesh(double radius, uint32_t subdivisionLevel = 0);
    
    // Clear the mesh (remove all vertices, edges and faces)
    void clear();
    
    // Create a base icosahedron
    void createIcosahedron(double radius);
    
    // Subdivide the current mesh
    void subdivide(uint32_t levels);
    
    // Validate the mesh structure and attempt to fix any inconsistencies
    void validateMesh();
    
    // Get the number of vertices in the mesh
    uint32_t getVertexCount() const { return static_cast<uint32_t>(m_vertices.size()); }
    
    // Get the position of a vertex
    double3 getVertexPosition(uint32_t index) const;
    
    // Set the position of a vertex
    void setVertexPosition(uint32_t index, const double3& position);
    
    // Get the number of faces in the mesh
    uint32_t getFaceCount() const { return static_cast<uint32_t>(m_faces.size()); }
    
    // Get the vertices of a face
    std::vector<uint32_t> getFaceVertices(uint32_t faceIndex) const;
    
    // Get all neighbors of a face (faces that share an edge)
    std::vector<uint32_t> getFaceNeighbors(uint32_t faceIndex) const;
    
    // Find all edges connected to a vertex
    std::vector<uint32_t> findVertexEdges(uint32_t vertexIndex) const;
    
    // Calculate the area of a face
    double calculateFaceArea(uint32_t faceIndex) const;
    
    // Calculate the normal vector of a face
    double3 calculateFaceNormal(uint32_t faceIndex) const;
    
    // Add a triangular face to the mesh (returns the face index)
    uint32_t addFace(uint32_t v1, uint32_t v2, uint32_t v3);
    
    // Add a vertex to the mesh (returns the vertex index)
    uint32_t addVertex(const double3& position);
    
    // Constants
    static constexpr uint32_t INVALID_INDEX = static_cast<uint32_t>(-1);
    
    // Structure to store vertex data
    struct Vertex {
        double3 position;    // Position in 3D space
        uint32_t edgeIndex;  // Index of one edge that starts from this vertex
        
        Vertex() : position(0, 0, 0), edgeIndex(INVALID_INDEX) {}
        Vertex(const double3& pos) : position(pos), edgeIndex(INVALID_INDEX) {}
    };
    
    // Structure to store edge data in the Winged-Edge format
    struct Edge {
        uint32_t startVertex; // Index of start vertex
        uint32_t endVertex;   // Index of end vertex
        uint32_t leftFace;    // Index of face on the left
        uint32_t rightFace;   // Index of face on the right
        
        // Edges around faces (classic Winged-Edge structure)
        uint32_t leftCW;      // Edge clockwise from this edge around left face
        uint32_t leftCCW;     // Edge counterclockwise from this edge around left face
        uint32_t rightCW;     // Edge clockwise from this edge around right face
        uint32_t rightCCW;    // Edge counterclockwise from this edge around right face
        
        Edge();
    };
    
    // Structure to store face data
    struct Face {
        uint32_t edgeIndex;   // Index of one edge that borders this face
        
        Face() : edgeIndex(INVALID_INDEX) {}
    };

private:
    // Private helper methods
    std::string edgeKey(uint32_t v1, uint32_t v2) const;
    Edge createEdge(uint32_t startVertex, uint32_t endVertex);
    uint32_t addEdge(uint32_t startVertex, uint32_t endVertex);
    uint32_t findEdge(uint32_t startVertex, uint32_t endVertex) const;
    uint32_t findOrCreateEdge(uint32_t startVertex, uint32_t endVertex);
    void linkEdges(uint32_t e1, uint32_t e2, bool startToStart, bool preserveWinding);
    void attachEdgeToFace(uint32_t edgeIndex, uint32_t faceIndex, bool isLeftFace);
    void completeEdgeLoop(uint32_t firstEdge, uint32_t lastEdge, uint32_t faceIndex);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2, std::unordered_map<std::string, uint32_t>& midpointCache, double sphereRadius);

    // Member variables
    std::vector<Vertex> m_vertices;  // Vertices of the mesh
    std::vector<Edge> m_edges;       // Edges of the mesh (Winged-Edge structure)
    std::vector<Face> m_faces;       // Triangular faces of the mesh
    
    // Edge lookup for faster edge finding during construction
    std::unordered_map<std::string, uint32_t> m_edgeMap;
};

