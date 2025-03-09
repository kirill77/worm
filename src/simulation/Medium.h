#pragma once

#include <memory>
#include <vector>
#include <string>
#include <array>
#include "Protein.h"
#include "MRNA.h"
#include "math/vector.h"

class Medium
{
private:
    // Represents contents of one grid cell
    struct GridCell {
        std::vector<ProteinPopulation> m_proteins;
        std::vector<std::shared_ptr<MRNA>> m_pMRNAs;
        double m_fAtp;  // ATP level in this grid cell
        
        // Helper to find or create protein population
        ProteinPopulation& findOrCreatePopulation(const ProteinPopulation& protein);
        
        GridCell() : m_fAtp(0.0) {}
    };

    static constexpr uint32_t GRID_RES = 3;  // 3x3x3 grid
    std::array<GridCell, GRID_RES * GRID_RES * GRID_RES> m_grid;
    
    static constexpr double DIFFUSION_RATE = 0.1;          // Rate of movement between cells
    static constexpr double CORTEX_BINDING_RATE = 0.3;     // Rate of binding to cortex
    static constexpr double CORTEX_UNBINDING_RATE = 0.1;   // Rate of unbinding from cortex
    static constexpr double ATP_DIFFUSION_RATE = 0.2;      // Rate of ATP diffusion between cells

public:
    static constexpr double MAX_ATP_PER_CELL = 100.0;      // Maximum ATP per grid cell

    // Add protein population to specific location
    void addProtein(const ProteinPopulation& protein, const float3& position);
    
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
    // Helper functions
    GridCell& findCell(const float3& position);
    const GridCell& findCell(const float3& position) const;
    bool isCortexCell(size_t index) const;  // Determines if a cell is on the boundary
    std::vector<size_t> getNeighborIndices(size_t cellIndex) const;
    
    // Update functions
    void updateProteinDiffusion(double dt);
    void updatePARDynamics(double dt);
    void translateMRNAs(double dt);
    void updateATPDiffusion(double dt);
    
    // Convert between grid indices and 3D coordinates
    uint32_t positionToIndex(const float3& position) const;
    float3 indexToPosition(size_t index) const;
    
    // Helper function to check protein type
    bool isPARProtein(const Protein& protein) const;
};

