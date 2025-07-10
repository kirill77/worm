#pragma once

#include <memory>
#include "Medium.h"
#include "Organelle.h"
#include "geometry/vectors/vector.h"
#include "physics/tensionSphere/tensionSphere.h"

/**
 * The Cortex class represents the cell membrane that separates
 * the internal cellular environment from the external environment.
 * It mediates interactions between internal and external media.
 */
struct Cortex : public Organelle
{
private:
    double m_fThickness; // Membrane thickness in micrometers
    TensionSphere m_tensionSphere;

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

    /**
     * Transport proteins from external to internal medium.
     * 
     * @param externalMedium Reference to the external medium
     * @param proteinName Name of the protein to transport
     * @param amount Amount of protein to transport
     * @param position Position where transport occurs
     * @return True if transport was successful
     */
    bool transportProteinInward(Medium& externalMedium, 
                               const std::string& proteinName, 
                               double amount, 
                               const float3& position);

    /**
     * Transport proteins from internal to external medium.
     * 
     * @param externalMedium Reference to the external medium
     * @param proteinName Name of the protein to transport
     * @param amount Amount of protein to transport
     * @param position Position where transport occurs
     * @return True if transport was successful
     */
    bool transportProteinOutward(Medium& externalMedium,
                                const std::string& proteinName,
                                double amount,
                                const float3& position);

    /**
     * Transport ATP from external to internal medium.
     * 
     * @param externalMedium Reference to the external medium
     * @param amount Amount of ATP to transport
     * @param position Position where transport occurs
     * @return True if transport was successful
     */
    bool transportATPInward(Medium& externalMedium,
                           double amount,
                           const float3& position);

    /**
     * Transport ATP from internal to external medium.
     * 
     * @param externalMedium Reference to the external medium
     * @param amount Amount of ATP to transport
     * @param position Position where transport occurs
     * @return True if transport was successful
     */
    bool transportATPOutward(Medium& externalMedium,
                            double amount,
                            const float3& position);

    // Getters and setters
    double getThickness() const { return m_fThickness; }
    void setThickness(double fThickness) { m_fThickness = fThickness; }

    const TensionSphere &getTensionSphere() const
    {
        return m_tensionSphere;
    }
}; 