#pragma once

#include <vector>
#include <memory>
#include <cstdint>
#include <array>
#include "../math/vector.h"
#include <unordered_map>
#include <string>

/**
 * @brief Represents a single cell in the geodesic sphere
 */
class SphereCell {
public:
    SphereCell();
    
    double getArea() const { return m_fArea; }
    void setArea(double fArea) { m_fArea = fArea; }
    
    double getAreaScaler() const { return m_fAreaScaler; }
    void setAreaScaler(double fScaler) { m_fAreaScaler = fScaler; }
    
    double getTension() const { return m_fTension; }
    void setTension(double fTension) { m_fTension = fTension; }
    
    double getTensionScaler() const { return m_fTensionScaler; }
    void setTensionScaler(double fScaler) { m_fTensionScaler = fScaler; }
    
    const std::vector<uint32_t>& getNeighbors() const { return m_neighbors; }
    void clearNeighbors() { m_neighbors.clear(); }
    void addNeighbor(uint32_t neighborIndex);

private:
    double m_fArea;           // Current area of the cell
    double m_fAreaScaler;     // Scaling factor for area (1.0 = normal)
    double m_fTension;        // Current tension in the cell
    double m_fTensionScaler;  // Scaling factor for tension (1.0 = normal)
    std::vector<uint32_t> m_neighbors; // Indices of neighboring cells
};

/**
 * @brief Class for simulating tension forces in a spherical cell cortex using a geodesic sphere model
 * implemented with Winged-Edge Data Structure for efficient mesh traversal
 */
class TensionSphere {
public:
    /**
     * @brief Constructor
     * @param subdivisionLevel Number of times to subdivide the base icosahedron (increases detail)
     */
    TensionSphere(uint32_t subdivisionLevel = 2);
    
    /**
     * @brief Get the total number of cells in the geodesic sphere
     */
    uint32_t getCellCount() const;
    
    /**
     * @brief Get a specific cell by index
     * @param index Index of the cell to retrieve
     * @return Reference to the cell
     */
    SphereCell& getCell(uint32_t index);
    
    /**
     * @brief Get a specific cell by index (const version)
     * @param index Index of the cell to retrieve
     * @return Const reference to the cell
     */
    const SphereCell& getCell(uint32_t index) const;
    
    /**
     * @brief Advance the simulation by one time step
     * @param fDtSec The time step size in seconds
     */
    void makeTimeStep(double fDtSec);
    
    /**
     * @brief Calculate and return the total tension energy in the system
     */
    double getTotalTensionEnergy() const;
    
    /**
     * @brief Reset all cells to balanced state (equal area and tension)
     */
    void resetToBalancedState();
    
    /**
     * @brief Set the stiffness coefficient for vertex movement
     * @param stiffness Stiffness value (higher = stiffer, less movement)
     */
    void setStiffness(double stiffness) { m_stiffness = stiffness; }
    
    /**
     * @brief Set the damping coefficient for vertex movement
     * @param damping Damping value (higher = more damping, less oscillation)
     */
    void setDamping(double damping) { m_damping = damping; }
    
    /**
     * @brief Get vertex position
     * @param index Vertex index
     * @return Position vector
     */
    double3 getVertexPosition(uint32_t index) const;
    
    /**
     * @brief Get the vertex count
     */
    uint32_t getVertexCount() const { return static_cast<uint32_t>(m_vertices.size()); }

private:
    // Forward declarations for Winged-Edge data structure components
    struct Vertex;
    struct Edge;
    struct Face;
    
    // Struct for representing a vertex in 3D space with velocity and force
    struct Vertex {
        double3 position;     // Position vector
        double3 velocity;     // Velocity vector
        double3 force;        // Force accumulator
        uint32_t edgeIndex;   // Index of one edge that starts from this vertex
    };
    
    // Struct for representing an edge in the Winged-Edge structure
    struct Edge {
        uint32_t startVertex; // Index of start vertex
        uint32_t endVertex;   // Index of end vertex
        uint32_t leftFace;    // Index of face on the left
        uint32_t rightFace;   // Index of face on the right
        
        // Edges around vertices and faces (clockwise/counterclockwise traversal)
        uint32_t startCW;     // Edge clockwise from this edge around start vertex
        uint32_t startCCW;    // Edge counterclockwise from this edge around start vertex
        uint32_t endCW;       // Edge clockwise from this edge around end vertex
        uint32_t endCCW;      // Edge counterclockwise from this edge around end vertex
        uint32_t leftCW;      // Edge clockwise from this edge around left face
        uint32_t leftCCW;     // Edge counterclockwise from this edge around left face
        uint32_t rightCW;     // Edge clockwise from this edge around right face
        uint32_t rightCCW;    // Edge counterclockwise from this edge around right face
    };
    
    // Struct for representing a triangular face
    struct Face {
        uint32_t edgeIndex;   // Index of one edge that borders this face
        double area;          // Area of the face
        double restArea;      // Rest area (equilibrium)
    };

    // Private methods for mesh construction and manipulation
    void createIcosahedron();
    void subdivide(uint32_t level);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2);
    void setupCellNeighbors();
    void normalizeVertex(Vertex& v);
    double calculateFaceArea(uint32_t faceIndex) const;
    void calculateForces();
    void integrateMotion(double fDtSec);
    void updateCellAreas();
    void enforceSphericalConstraint();
    double3 calculateFaceNormal(uint32_t faceIndex) const;
    
    // Winged-Edge specific methods
    Edge createEdge(uint32_t startVertex, uint32_t endVertex);
    void linkEdges(uint32_t e1, uint32_t e2, bool startToStart, bool preserveWinding);
    void attachEdgeToFace(uint32_t edgeIndex, uint32_t faceIndex, bool isLeftFace);
    void completeEdgeLoop(uint32_t firstEdge, uint32_t lastEdge, uint32_t faceIndex);
    std::vector<uint32_t> getFaceVertices(uint32_t faceIndex) const;
    uint32_t addEdge(uint32_t startVertex, uint32_t endVertex);
    uint32_t findEdge(uint32_t startVertex, uint32_t endVertex) const;
    uint32_t findOrCreateEdge(uint32_t startVertex, uint32_t endVertex);
    uint32_t addFace(uint32_t v1, uint32_t v2, uint32_t v3);
    
    // Helper to generate a unique string key for an edge
    std::string edgeKey(uint32_t v1, uint32_t v2) const;

    // Member variables
    std::vector<Vertex> m_vertices;  // Vertices of the geodesic sphere
    std::vector<Edge> m_edges;       // Edges of the geodesic sphere (Winged-Edge structure)
    std::vector<Face> m_faces;       // Triangular faces of the geodesic sphere
    std::vector<SphereCell> m_cells; // Cells (corresponding to faces)
    double m_simulationTime;         // Current simulation time in seconds
    
    // Edge lookup for faster edge finding during construction
    std::unordered_map<std::string, uint32_t> m_edgeMap;
    
    // Physics parameters
    double m_stiffness;       // Stiffness coefficient (resistance to deformation)
    double m_damping;         // Damping coefficient (reduces oscillations)
    double m_sphereRadius;    // Radius constraint for the sphere
    
    // Constants used in Winged-Edge structure
    static constexpr uint32_t INVALID_INDEX = static_cast<uint32_t>(-1);
};

