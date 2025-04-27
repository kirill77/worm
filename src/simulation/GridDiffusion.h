#pragma once

#include "Grid.h"
#include <vector>
#include <string>

class GridDiffusion
{
public:
    GridDiffusion();

    // Update diffusion for a specific molecule type
    void updateDiffusion(Grid& grid, double dt);

private:

    static constexpr double DIFFUSION_RATE = 0.1; // Rate of movement between cells

    // this function copies GridCell::m_fAtp into GridCell::m_proteins
    void copyATPToProteins(Grid& grid) const;
    // this function copies GridCell::m_proteins into GridCell::m_fAtp
    void copyATPFromProteins(Grid& grid) const;

    // Helper function to compute diffusion amount
    double computeDiffusionAmount(double moleculeCount, size_t numNeighbors, double dt) const;
}; 