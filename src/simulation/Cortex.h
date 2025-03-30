#pragma once

#include <memory>
#include "Medium.h"
#include "ProteinBindingSurface.h"
#include "math/vector.h"

/**
 * The Membrane class represents the cell membrane that separates
 * the internal cellular environment from the external environment.
 * It mediates interactions between internal and external media.
 */
class Cortex : public ProteinBindingSurface
{
private:
    std::shared_ptr<Medium> m_pInternalMedium;  // Internal cellular environment
    double m_fThickness;                        // Membrane thickness in micrometers

public:
    /**
     * Constructor that initializes the membrane with an internal medium.
     * The external medium is not stored and must be passed to methods that need it.
     * 
     * @param pInternalMedium Shared pointer to the cell's internal medium
     * @param fThickness Membrane thickness in micrometers
     * @param fSurfaceArea Surface area in square micrometers
     */
    Cortex(std::shared_ptr<Medium> pInternalMedium,
        double fThickness = 0.01); // Default 10nm thickness

    /**
     * Update the membrane state.
     * This method updates membrane dynamics and can be extended to include
     * passive transport, signal transduction, etc.
     * 
     * @param dt Time step in seconds
     */
    void update(double dt);
    
    /**
     * Initialize binding sites in the internal medium.
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
    std::shared_ptr<Medium> getInternalMedium() const { return m_pInternalMedium; }
    Medium& getInternalMedium() { return *m_pInternalMedium; }
    
    double getThickness() const { return m_fThickness; }
    void setThickness(double fThickness) { m_fThickness = fThickness; }
}; 