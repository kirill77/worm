#pragma once

#include <vector>
#include <memory>
#include <assert.h>
#include "CellTypes.h"
#include "Chromosome.h"
#include "chemistry/molecules/StringDict.h"
#include "chemistry/molecules/Molecule.h"

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
    static constexpr double fMRNA_EXPORT = 2.0;            // Cost per mRNA export from nucleus
};

class BVHMesh;

class Cell : public std::enable_shared_from_this<Cell>
{
private:
    std::vector<std::shared_ptr<class Organelle>> m_pOrganelles;
    std::shared_ptr<class Medium> m_pInternalMedium;  // Internal cellular environment
    CellCycleState m_cellCycleState;
    CellType m_type;  // Store type just for spindle creation
    std::vector<Chromosome> m_chromosomes;  // Store chromosomes for delayed organelle creation
    Species m_species = Species::GENERIC;   // Biological species of the cell

    // Helper functions
    void checkCellCycleTransitions();
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
    Cell(std::shared_ptr<Medium> pInternalMedium, const std::vector<Chromosome>& chromosomes, CellType type = CellType::Zygote, Species species = Species::GENERIC);

public:
    // Static factory method to create a cell
    static std::shared_ptr<Cell> createCell(std::shared_ptr<Medium> pInternalMedium,
                                          const std::vector<Chromosome>& chromosomes, 
                                          CellType type = CellType::Zygote,
                                          Species species = Species::GENERIC);
    
    void update(double fDt);
    CellCycleState getCellCycleState() const { return m_cellCycleState; }

    Species getSpecies() const { return m_species; }
    
    // Access to internal medium
    Medium& getInternalMedium() const { return *m_pInternalMedium; }
    
    // Organelle management
    void addOrganelle(StringDict::ID id, std::shared_ptr<Organelle> pOrganelle);
    std::shared_ptr<Organelle> getOrganelle(StringDict::ID id) const;  // Generic organelle getter

    // ATP-related functions
    bool consumeATP(double fAmount);
    
    // Cortex BVH moved to Cortex; queries should go through Cortex
};

