#pragma once

#include <memory>
#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include "chemistry/molecules/Molecule.h"
#include "biology/organelles/CortexLocation.h"

#include "geometry/vectors/vector.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "chemistry/molecules/GridCell.h"
#include "Grid.h"
#include "GridDiffusion.h"
#include "chemistry/interactions/ResourceDistributor.h"

// Forward declaration to avoid circular include
class Cortex;

class Medium
{
private:
    Grid m_grid;
    GridDiffusion m_diffusion;
    double m_fVolumeMicroM;  // Volume in micrometers

    static constexpr double ATP_DIFFUSION_RATE = 0.2;      // Rate of ATP diffusion between cells
    static constexpr int DIFFUSION_SAMPLES = 1000;         // Number of random samples per diffusion update
    static constexpr double DIFFUSION_SIGMA = 0.2;         // Standard deviation for diffusion distance (as fraction of medium size)

public:
    static constexpr double MAX_ATP_PER_CELL = 1e10;      // Maximum ATP per grid cell

    Medium();

    // Add molecule population to specific location
    void addMolecule(const MPopulation& population, const float3& position);
    // Move molecules from grid cells into the provided binding sites based on m_normalized
    // Only molecules listed in bindableMolecules are transferred
    void toBindingSites(std::vector<CortexMolecules>& bindingSites, const std::vector<Molecule>& bindableMolecules);
    
    // ATP-related methods
    void addATP(double amount, const float3& position);
    bool consumeATP(double amount, const float3& position);
    double getAvailableATP(const float3& position) const;
    
    // Get molecule number at a specific location
    double getMoleculeNumber(const Molecule& molecule, const float3& position) const;
    // Get molecule concentration (per µm^3) at a specific location
    double getMoleculeConcentration(const Molecule& molecule, const float3& position) const;
    
    
    // Get volume in micrometers
    double getVolumeMicroM() const { return m_fVolumeMicroM; }

    // Update per-cell volumes based on cortex world mapping
    void updateGridCellVolumes(Cortex& cortex);
    
    // Main update function
    void update(double dt);

private:
    ResourceDistributor m_resDistributor;

    // Update functions
    void updateMoleculeInteraction(double dt);
};

