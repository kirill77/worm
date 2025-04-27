#pragma once

#include <memory>
#include <vector>
#include <string>
#include <array>
#include <unordered_map>
#include "Molecule.h"
#include "MRNA.h"
#include "math/vector.h"
#include "ProteinWiki.h"
#include "GridCell.h"
#include "Grid.h"
#include "ResourceDistributor.h"

class Medium
{
private:
    Grid m_grid;

    static constexpr double DIFFUSION_RATE = 0.1;          // Rate of movement between cells
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
    
    // ATP-related methods
    void addATP(double amount, const float3& position);
    bool consumeATP(double amount, const float3& position);
    double getAvailableATP(const float3& position) const;
    
    // Get number of proteins at a specific location
    double getProteinNumber(const std::string& proteinName, const float3& position) const;
    
    // Get total number of proteins across all cells
    double getTotalProteinNumber(const std::string& proteinName) const;
    
    // Main update function
    void update(double dt);

private:
    ResourceDistributor m_resDistributor;

    // Update functions
    void updateProteinDiffusion(double dt);
    void updateProteinInteraction(double dt);
    void translateMRNAs(double dt);
    void updateATPDiffusion(double dt);
    
    // Helper functions for physically-based diffusion
    float3 generateRandomPosition() const;
    float3 generateRandomDirection() const;
    float generateRandomDistance(double dt) const;
};

