#pragma once

#include <memory>
#include "Medium.h"
#include "Organelle.h"
#include "geometry/vectors/vector.h"

/**
 * The Cortex class represents the cell membrane that separates
 * the internal cellular environment from the external environment.
 * It mediates interactions between internal and external media.
 */
struct Cortex : public Organelle
{
private:
    double m_fThickness; // Membrane thickness in micrometers
    std::shared_ptr<class BVHMesh> m_pCortexBVH;
    std::shared_ptr<class TensionSphere> m_pTensionSphere;

public:
    /**
     * Constructor that initializes the cortex.
     * 
     * @param pCell Weak pointer to the cell containing this cortex
     * @param fThickness Membrane thickness in micrometers
     */
    Cortex(std::weak_ptr<Cell> pCell, double fThickness = 0.01); // Default 10nm thickness

    /**
     * Update the cortex state.
     * This method updates cortex dynamics and can be extended to include
     * passive transport, signal transduction, etc.
     * 
     * @param fDtSec Time step in seconds
     * @param cell Reference to the cell containing this cortex
     */
    void update(double fDtSec, Cell& cell) override;
    
    /**
     * Initialize binding sites in the cell's internal medium.
     * This creates binding sites throughout the medium
     * that allow proteins to bind to the cell membrane surface.
     * 
     * @param totalAmount Total amount of binding sites to distribute
     * @return True if binding sites were successfully added
     */
    bool initializeBindingSites(double totalAmount = 1000.0);

    // Getters and setters
    double getThickness() const { return m_fThickness; }
    void setThickness(double fThickness) { m_fThickness = fThickness; }

    // No external accessors for simulation-owned resources; managed internally in update()
}; 