#pragma once
#include <memory>
#include "molecules/BindingSurface.h"

class Cell;
class Medium;
enum class CellCycleState;

/**
 * Base class for all cellular organelles.
 * Organelles can bind proteins and interact with the cell's internal medium.
 */
class Organelle : public BindingSurface
{
private:
    std::weak_ptr<Cell> m_pCell;  // Reference to the cell containing this organelle

public:
    /**
     * Constructor that initializes the organelle with a reference to its cell.
     * 
     * @param pCell Weak pointer to the cell containing this organelle
     */
    Organelle(std::weak_ptr<Cell> pCell)
        : BindingSurface()
        , m_pCell(pCell)
    {}
    
    /**
     * Base update function that takes cell
     * 
     * @param dt Time step in seconds
     * @param cell Reference to the cell containing this organelle
     */
    virtual void update(double dt, Cell& cell) = 0;
    
    /**
     * Virtual destructor for proper cleanup
     */
    virtual ~Organelle() {}

    /**
     * Get the cell containing this organelle
     * 
     * @return Shared pointer to the cell, or nullptr if the cell no longer exists
     */
    std::shared_ptr<Cell> getCell() const { return m_pCell.lock(); }
};

