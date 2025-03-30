#pragma once

#include <vector>
#include <memory>
#include <cstdint>
#include <array>
#include "math/vector.h"

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
    // Struct for representing a vertex in 3D space with velocity and force
    struct Vertex {
        double3 position;    // Position vector
        double3 velocity;    // Velocity vector
        double3 force;       // Force accumulator
    };
    
    // Struct for representing a triangular face
    struct Face {
        uint32_t v1, v2, v3;  // Vertex indices
        double area;          // Area of the face
        double restArea;      // Rest area (equilibrium)
    };

    // Private methods
    void createIcosahedron();
    void subdivide(uint32_t level);
    uint32_t getMidpoint(uint32_t v1, uint32_t v2);
    void setupCellNeighbors();
    void normalizeVertex(Vertex& v);
    double calculateFaceArea(const Face& face) const;
    void calculateForces();
    void integrateMotion(double fDtSec);
    void updateCellAreas();
    void enforceSphericalConstraint();
    double3 calculateFaceNormal(const Face& face) const;

    // Member variables
    std::vector<Vertex> m_vertices;  // Vertices of the geodesic sphere
    std::vector<Face> m_faces;       // Triangular faces of the geodesic sphere
    std::vector<SphereCell> m_cells; // Cells (corresponding to faces)
    double m_simulationTime;         // Current simulation time in seconds
    
    // Physics parameters
    double m_stiffness;       // Stiffness coefficient (resistance to deformation)
    double m_damping;         // Damping coefficient (reduces oscillations)
    double m_sphereRadius;    // Radius constraint for the sphere
};

