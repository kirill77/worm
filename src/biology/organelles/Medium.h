#pragma once

#include <memory>
#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include "chemistry/Molecule.h"
#include "chemistry/MRNA.h"
#include "chemistry/TRNA.h"
#include "geometry/vectors/vector.h"
#include "chemistry/MoleculeWiki.h"
#include "chemistry/GridCell.h"
#include "Grid.h"
#include "GridDiffusion.h"
#include "chemistry/ResourceDistributor.h"

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

    // Add protein population to specific location
    void addProtein(const MPopulation& protein, const float3& position);
    
    // Add mRNA to specific location
    void addMRNA(std::shared_ptr<MRNA> mRNA, const float3& position);
    
    // Add tRNA to specific location
    void addTRNA(std::shared_ptr<TRNA> tRNA, const float3& position);
    
    // ATP-related methods
    void addATP(double amount, const float3& position);
    bool consumeATP(double amount, const float3& position);
    double getAvailableATP(const float3& position) const;
    
    // Get number of proteins at a specific location
    double getProteinNumber(const std::string& proteinName, const float3& position) const;
    
    // Get total number of proteins across all cells
    double getTotalProteinNumber(const std::string& proteinName) const;
    
    // Get volume in micrometers
    double getVolumeMicroM() const { return m_fVolumeMicroM; }
    
    // Main update function
    void update(double dt);

private:
    ResourceDistributor m_resDistributor;

    // Update functions
    void updateProteinInteraction(double dt);
    void translateMRNAs(double dt);
};

