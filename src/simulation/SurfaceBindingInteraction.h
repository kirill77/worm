#pragma once

#include "ProteinInteraction.h"
#include "GridCell.h"

/**
 * Represents the binding of proteins to a protein binding surface.
 * This interaction allows proteins to associate with any surface that
 * has injected binding sites.
 */
class SurfaceBindingInteraction : public ProteinInteraction
{
public:
    struct Parameters {
        double bindingRate;         // Rate at which proteins bind to surface sites
        double dissociationRate;    // Rate at which proteins dissociate from surface sites
        double saturationConstant;  // For binding kinetics (used in denominator)
    };
    
    // Constructor
    SurfaceBindingInteraction(
        const std::string& proteinName, 
        const std::string& bindingSiteName,
        const std::string& boundComplexName,
        const Parameters& params);
    
    // Apply binding to proteins in the cell
    bool apply(GridCell& cell, double dt, double& atpConsumed) const override;
    
private:
    std::string m_proteinName;            // Name of the protein to bind
    std::string m_bindingSiteName;        // Name of binding site proteins
    std::string m_boundComplexName;       // Name of protein bound to surface
    double m_bindingRate;
    double m_dissociationRate;
    double m_saturationConstant;
}; 