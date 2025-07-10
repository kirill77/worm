#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "geometry/vectors/vector.h"
#include "geometry/edgeMesh/edgeMesh.h"

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

    void makeTimeStep(double fDtSec);

    /**
     * @brief Get the underlying EdgeMesh
     * @return Shared pointer to the EdgeMesh
     */
    std::shared_ptr<EdgeMesh> getEdgeMesh() const { return m_pMesh; }

    /**
     * @brief Get the volume of the tension sphere
     * @return Volume value
     */
    double getVolume() const;

    /**
     * @brief Set the volume of the tension sphere
     * @param volume Volume value to set
     */
    void setVolume(double volume);

private:
    // The underlying mesh data structure
    std::shared_ptr<EdgeMesh> m_pMesh;

    // Velocity for each vertex (same indexing as mesh vertices)
    std::vector<double3> m_vertexVelocities;

    // Edge rest lengths (computed from initial mesh)
    std::vector<double> m_edgeRestLengths;
    
    // Edge connectivity (pairs of vertex indices)
    std::vector<std::pair<uint32_t, uint32_t>> m_edgeConnectivity;

    // constants controlling spring behaviour
    double m_fSpringC = 0.1, m_fDampingCoeff = 1;

    // Volume of the tension sphere
    double m_fVolume;

    // Helper methods
    void initializePhysics();
    void computeSpringForces(std::vector<double3>& forces, double dt);
    void integrateMotion(const std::vector<double3>& forces, double dt);
};

