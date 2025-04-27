#pragma once

#include "Grid.h"
#include <vector>
#include <string>

class GridDiffusion
{
public:
    struct Parameters {
        double diffusionRate;      // Rate of movement between cells
        int diffusionSamples;      // Number of random samples per diffusion update
        double diffusionSigma;     // Standard deviation for diffusion distance
    };

    GridDiffusion(const Parameters& params);

    // Update diffusion for a specific molecule type
    void updateDiffusion(Grid& grid, const std::string& moleculeName, double dt);

private:
    // this function copies GridCell::m_fAtp into GridCell::m_proteins
    void copyATPToProteins(Grid& grid) const;
    // this function copies GridCell::m_proteins into GridCell::m_fAtp
    void copyATPFromProteins(Grid& grid) const;

    Parameters m_params;

    // Helper function to compute diffusion amount
    double computeDiffusionAmount(double moleculeCount, size_t numNeighbors, double dt) const;
}; 