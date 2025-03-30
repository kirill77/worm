#pragma once

#include <vector>
#include <memory>
#include "Cortex.h"
#include "CellTypes.h"
#include "Chromosome.h"

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
    static constexpr double fPROTEIN_SYNTHESIS = 4.0;      // Cost per protein molecule
    static constexpr double fCHROMOSOME_CONDENSATION = 10.0; // Cost per chromosome
    static constexpr double fSPINDLE_FORMATION = 15.0;     // Cost for mitotic spindle
    static constexpr double fCHROMOSOME_MOVEMENT = 5.0;    // Cost per chromosome per second during anaphase
    static constexpr double fMEMBRANE_FUSION = 8.0;        // Cost for membrane fusion events
    static constexpr double fMRNA_SYNTHESIS = 2.0;         // Cost per mRNA molecule
};

class Cell
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<Cortex> m_pCortex;
    CellCycleState m_cellCycleState;
    CellType m_type;  // Store type just for spindle creation

    // Helper functions
    void checkCellCycleTransitions();
    std::shared_ptr<class Mitochondrion> getMitochondrion() const;
    void createSpindle();
    void destroySpindle();

public:
    // Constructor that takes a membrane instead of a medium
    Cell(std::shared_ptr<Cortex> pCortex, const std::vector<Chromosome>& chromosomes, CellType type = CellType::Zygote);
    
    void update(double fDt);
    CellCycleState getCellCycleState() const { return m_cellCycleState; }
    
    // Access to internal medium
    Medium& getInternalMedium() const { return m_pCortex->getInternalMedium(); }
    
    std::shared_ptr<class Spindle> getSpindle() const;  // Made public for Chromosome access

    // ATP-related functions
    bool consumeATP(double fAmount);
};

