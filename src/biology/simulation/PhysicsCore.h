#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/edgeMesh.h"
#include "physics/ForceGenerator.h"
#include "physics/PhysicsConstraints.h"
#include "physics/PhysicsIntegrator.h"
#include "SoftBodyMeshAdapter.h"

class Cell;

/**
 * Core physics simulator for the cell that advances mechanics and constraints
 * over shared mesh representations (e.g., cortex mesh, microtubules via adapters).
 */
class PhysicsCore {
public:
    PhysicsCore() = default;

    // Initialize physics core with cell reference; pulls mesh and volume from cell
    void initialize(std::shared_ptr<Cell> pCell);

    void makeTimeStep(double fDtSec);

private:
    // Reference to the cell (for accessing cortex mesh and medium volume)
    std::shared_ptr<Cell> m_pCell;

    // Underlying mesh (currently cortex); more adapters can be added later
    std::shared_ptr<EdgeMesh> m_pCortexMesh;

    // Edge rest lengths (computed from initial mesh)
    std::vector<double> m_edgeRestLengths;

    // Mesh adapter (reused across timesteps to avoid repeated allocations)
    std::shared_ptr<SoftBodyMeshAdapter> m_pMeshAdapter;

    // constants controlling spring behaviour
    double m_fSpringC = 0.1, m_fDampingCoeff = 1;

    // Target volume
    double m_fVolume;

    // Physics integrator managing all bodies
    PhysicsIntegrator m_integrator;

    // Pluggable force generators acting on the mesh
    std::vector<std::unique_ptr<IForceGenerator>> m_forceGenerators;

    // Per-body constraints (XPBD-style)
    std::vector<std::unique_ptr<IConstraint>> m_constraints;

    // Internal initialization helpers
    void initializePhysics();
};


