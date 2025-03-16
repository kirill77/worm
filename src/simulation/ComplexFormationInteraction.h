#pragma once

#include "ProteinInteraction.h"
#include "GridCell.h"

/**
 * Represents a complex formation interaction where two proteins
 * bind together to form a functional complex.
 */
class ComplexFormationInteraction : public ProteinInteraction
{
public:
    struct Parameters {
        double bindingRate;        // Rate at which proteins form complexes
        double dissociationRate;   // Rate at which complexes break apart
        double saturationConstant; // Saturation constant for binding kinetics
        std::string complexName;   // Name of the resulting complex
    };
    
    // Constructor
    ComplexFormationInteraction(const std::string& proteinA, 
                              const std::string& proteinB, 
                              const Parameters& params);
    
    // Apply complex formation to proteins in the cell
    bool apply(GridCell& cell, double dt, double& atpConsumed) const override;

private:
    double m_bindingRate;
    double m_dissociationRate;
    double m_saturationConstant;
    std::string m_complexName;
}; 