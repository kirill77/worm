#pragma once

#include <vector>
#include <memory>
#include "Medium.h"

enum class CellCycleState
{
    INTERPHASE,
    PROPHASE,
    METAPHASE,
    ANAPHASE,
    TELOPHASE,
    CYTOKINESIS
};

// ATP costs for various cellular processes
struct ATPCosts
{
    static constexpr double PROTEIN_SYNTHESIS = 4.0;      // Cost per protein molecule
    static constexpr double CHROMOSOME_CONDENSATION = 10.0; // Cost per chromosome
    static constexpr double SPINDLE_FORMATION = 15.0;     // Cost for mitotic spindle
    static constexpr double CHROMOSOME_MOVEMENT = 5.0;    // Cost per chromosome per second during anaphase
    static constexpr double MEMBRANE_FUSION = 8.0;        // Cost for membrane fusion events
    static constexpr double MRNA_SYNTHESIS = 2.0;         // Cost per mRNA molecule
};

class Cell
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<Medium> m_pMedium;
    CellCycleState m_cellCycleState;

    // Helper function to check conditions for cell cycle transitions
    void checkCellCycleTransitions();
    std::shared_ptr<class Mitochondrion> getMitochondrion() const;

public:
    Cell(std::shared_ptr<Medium> pMedium);
    void update(double dt);
    CellCycleState getCellCycleState() const { return m_cellCycleState; }
    std::shared_ptr<Medium> getMedium() const { return m_pMedium; }

    // ATP-related functions
    double getAvailableATP() const;
    bool consumeATP(double amount);  // Returns false if not enough ATP
};

