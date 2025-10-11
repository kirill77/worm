#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/edgeMesh.h"
#include "ForceGenerator.h"
#include "PhysicsConstraints.h"

/**
 * @brief Class for simulating tension forces in a spherical cell cortex using a geodesic sphere model
 * implemented with Winged-Edge Data Structure for efficient mesh traversal
 */
class TensionSphere {
public:
    /**
     * @brief Constructor
     * @param subdivisionLevel Number of times to subdivide the base icosahedron (increases detail)
     * @param volume Initial volume of the tension sphere
     */
    TensionSphere(uint32_t subdivisionLevel = 2, double volume = 0.0);

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

    // constants controlling spring behaviour
    double m_fSpringC = 0.1, m_fDampingCoeff = 1;

    // Volume of the tension sphere
    double m_fVolume;

    // Pluggable force generators acting on the mesh
    std::vector<std::unique_ptr<IForceGenerator>> m_forceGenerators;

    // Per-body constraints (XPBD-style)
    std::vector<std::unique_ptr<IConstraint>> m_constraints;

    // Helper methods
    void initializePhysics();
};

