#pragma once
#include <memory>
#include "ProteinBindingSurface.h"

class Cell;
class Medium;
enum class CellCycleState;

/**
 * Base class for all cellular organelles.
 * Organelles can bind proteins and interact with the cell's internal medium.
 */
class Organelle : public ProteinBindingSurface
{
public:
    /**
     * Constructor that initializes the organelle with a specific surface area.
     * 
     * @param fSurfaceArea Surface area in square micrometers
     */
    Organelle(double fSurfaceArea = 1.0)
        : ProteinBindingSurface(fSurfaceArea)  // Default binding capacity of 1e6 proteins per square micrometer
    {}
    
    /**
     * Base update function that takes cell and medium
     * 
     * @param dt Time step in seconds
     * @param cell Reference to the cell containing this organelle
     * @param medium Reference to the cell's internal medium
     */
    virtual void update(double dt, Cell& cell, Medium& medium) = 0;
    
    /**
     * Virtual destructor for proper cleanup
     */
    virtual ~Organelle() {}
};

