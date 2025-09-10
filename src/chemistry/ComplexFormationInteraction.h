#pragma once

#include "MoleculeInteraction.h"
#include "GridCell.h"
#include "Molecule.h"
#include "StringDict.h"

/**
 * Represents a complex formation interaction where two proteins
 * bind together to form a functional complex.
 */
class ComplexFormationInteraction : public MoleculeInteraction
{
public:
    struct Parameters {
        double bindingRate;        // Rate at which proteins form complexes
        double dissociationRate;   // Rate at which complexes break apart
        double saturationConstant; // Saturation constant for binding kinetics
        StringDict::ID complexId;  // ID of the resulting complex
    };
    
    // Constructor
    ComplexFormationInteraction(const Molecule& firstProtein, 
                              const Molecule& secondProtein, 
                              const Parameters& params);
    
    // Apply complex formation to proteins in the cell
    bool apply(GridCell& cell, double dt, ResourceDistributor& resDistributor) const override;
    
private:
    Molecule m_firstProtein;    // First protein in complex
    Molecule m_secondProtein;   // Second protein in complex
    double m_bindingRate;
    double m_dissociationRate;
    double m_saturationConstant;
    StringDict::ID m_complexId;
}; 