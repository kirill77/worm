#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "../math/vector.h"
#include "connectedMesh/connectedMesh.h"

/**
 * @brief Represents a single cell in the geodesic sphere
 */
class SphereCell {
public:
    SphereCell();
    
    double getArea() const { return m_fArea; }
    void setArea(double fArea) { m_fArea = fArea; }
    
    double getAreaScaler() const { return m_fAreaScaler; }
    void setAreaScaler(double fAreaScaler) { m_fAreaScaler = fAreaScaler; }
    
    double getTension() const { return m_fTension; }
    void setTension(double fTension) { m_fTension = fTension; }
    
    double getTensionScaler() const { return m_fTensionScaler; }
    void setTensionScaler(double fTensionScaler) { m_fTensionScaler = fTensionScaler; }
    
    const std::vector<uint32_t>& getNeighbors() const { return m_neighbors; }
    void clearNeighbors() { m_neighbors.clear(); }
    void addNeighbor(uint32_t neighborIndex);

private:
    double m_fArea;            // Current area of the cell
    double m_fAreaScaler;      // Multiplier for area calculations (default 1.0)
    double m_fTension;         // Current tension value
    double m_fTensionScaler;   // Multiplier for tension calculations (default 1.0)
    std::vector<uint32_t> m_neighbors;  // Indices of neighboring cells
};

/**
 * @brief Class for simulating tension forces in a spherical cell cortex using a geodesic sphere model
 * implemented with Winged-Edge Data Structure for efficient mesh traversal
 */
class TensionSphere {
public:
    /**
     * @brief Extended Vertex structure to include physics properties
     */
    struct Vertex {
        double3 velocity;  // Velocity vector of the vertex
        double3 force;     // Force accumulator for this vertex
        
        Vertex() : velocity(0.0), force(0.0) {}
    };
    
    /**
     * @brief Face extension to store rest area
     */
    struct Face {
        double area;      // Current area
        double restArea;  // Rest area (no tension when current == rest)
        
        Face() : area(0.0), restArea(0.0) {}
    };

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
    uint32_t getVertexCount() const { return m_pMesh->getVertexCount(); }
    
    /**
     * @brief Get the underlying ConnectedMesh
     * @return Shared pointer to the ConnectedMesh
     */
    std::shared_ptr<ConnectedMesh> getConnectedMesh() const { return m_pMesh; }

private:
    // The underlying mesh data structure
    std::shared_ptr<ConnectedMesh> m_pMesh;
    
    // Physics and simulation properties
    std::vector<Vertex> m_vertexData;       // Physics data for vertices
    std::vector<Face> m_faceData;           // Extended face data
    std::vector<SphereCell> m_cells;        // Cell data for each face
    
    double m_simulationTime;                // Current simulation time
    double m_stiffness;                     // Spring stiffness constant
    double m_damping;                       // Velocity damping factor (0-1)
    double m_sphereRadius;                  // Radius of the sphere
    
    // Helper functions
    void calculateForces();                 // Calculate forces on vertices based on tensions
    void integrateMotion(double fDtSec);    // Update positions based on forces
    void enforceSphericalConstraint();      // Ensure vertices stay on the sphere
    void updateCellAreas();                 // Update cell areas after vertex movement
    void setupCellNeighbors();              // Setup neighbor relationships between cells
    
    // Utility functions
    double3 calculateFaceNormal(uint32_t faceIndex) const;
};

