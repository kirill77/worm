#pragma once

#include "Organelle.h"
#include <memory>

class Cell;

class EReticulum : public Organelle
{
public:
    EReticulum(std::weak_ptr<Cell> pCell);
    void update(double dt, Cell& cell) override;

private:
    // Constants for ER behavior
    static constexpr double PROTEIN_SYNTHESIS_RATE = 1000.0;  // Proteins per second
    static constexpr double LIPID_SYNTHESIS_RATE = 500.0;     // Lipids per second
    static constexpr double ATP_COST_PER_PROTEIN = 4.0;       // ATP cost per protein synthesis
    static constexpr double ATP_COST_PER_LIPID = 2.0;         // ATP cost per lipid synthesis

    // Helper functions
    void synthesizeProteins(Medium& medium);
    void synthesizeLipids(Medium& medium);
    float3 generateRandomPosition() const;
}; 