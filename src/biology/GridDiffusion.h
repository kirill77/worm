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

    // Helper function to compute diffusion amount
    double computeDiffusionAmount(double moleculeCount, size_t numNeighbors, double dt) const;
}; 