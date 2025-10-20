#pragma once

#include <vector>
#include <array>
#include <unordered_map>
#include <string>
#include <memory>
#include "geometry/vectors/vector.h"
#include "geometry/mesh/edgeMesh.h"
#include "physics/PhysicsIntegrator.h"
#include "physics/VolumeConstraint.h"
#include "physics/PhysicsMesh.h"

class Cell;
class PhysCentrosome;

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

    // Cortex physics mesh (reused across timesteps to avoid repeated allocations)
    std::shared_ptr<PhysicsMesh> m_pCortexAdapter;

    // Centrosomes for dynein force calculations
    std::vector<std::shared_ptr<PhysCentrosome>> m_centrosomes;

    // constants controlling spring behaviour
    double m_fSpringC = 0.1, m_fDampingCoeff = 1;

    // Physics integrator managing complete physics pipeline
    PhysicsIntegrator m_integrator;

    // Volume constraint for dynamic volume updates
    std::shared_ptr<VolumeConstraintXPBD> m_pVolumeConstraint;
};


