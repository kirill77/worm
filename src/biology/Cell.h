#pragma once

#include <vector>
#include <memory>
#include <assert.h>
#include "CellTypes.h"
#include "Chromosome.h"
#include "chemistry/StringDict.h"

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

class Cell : public std::enable_shared_from_this<Cell>
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<class Medium> m_pInternalMedium;  // Internal cellular environment
    CellCycleState m_cellCycleState;
    CellType m_type;  // Store type just for spindle creation
    std::vector<Chromosome> m_chromosomes;  // Store chromosomes for delayed organelle creation

    // Helper functions
    void checkCellCycleTransitions();
    std::shared_ptr<class Mitochondrion> getMitochondrion() const;
    void createSpindle();
    void destroySpindle();
    void initializeOrganelles();  // Initialize organelles after construction
    void initializeCortex();      // Initialize cortex after construction
    
    // Organelle indexing helper
    size_t getOrganelleIndex(StringDict::ID id) const {
        assert(id >= StringDict::ID::ORGANELLE_START && id < StringDict::ID::ORGANELLE_END);
        return static_cast<size_t>(id) - static_cast<size_t>(StringDict::ID::ORGANELLE_START);
    }

    // Private constructor
    Cell(std::shared_ptr<Medium> pInternalMedium, const std::vector<Chromosome>& chromosomes, CellType type = CellType::Zygote);

public:
    // Static factory method to create a cell
    static std::shared_ptr<Cell> createCell(std::shared_ptr<Medium> pInternalMedium,
                                          const std::vector<Chromosome>& chromosomes, 
                                          CellType type = CellType::Zygote);
    
    void update(double fDt);
    CellCycleState getCellCycleState() const { return m_cellCycleState; }
    
    std::shared_ptr<class Cortex> getCortex() const;
    
    // Access to internal medium
    Medium& getInternalMedium() const { return *m_pInternalMedium; }
    
    std::shared_ptr<class Spindle> getSpindle() const;  // Made public for Chromosome access

    // Organelle management
    void addOrganelle(StringDict::ID id, std::shared_ptr<Organelle> pOrganelle);
    std::shared_ptr<class Centrosome> getCentrosome() const;  // Get centrosome if it exists

    // ATP-related functions
    bool consumeATP(double fAmount);
};

