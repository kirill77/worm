#pragma once

#include "ProteinInteraction.h"

/**
 * Represents a synergistic interaction between proteins
 * where proteins enhance each other's function or stability.
 */
class SynergisticInteraction : public ProteinInteraction
{
public:
    // Parameters specific to synergistic interactions
    struct Parameters {
        double enhancementFactor;        // How much one protein enhances the other
        double enhancementDuration;      // How long the enhancement lasts
        double saturationConstant;       // Saturation constant for Hill-type kinetics
        ProteinInteraction::Mechanism mechanism;   // Mechanism of synergy
        double atpCost;                  // ATP consumed per unit of enhancement
    };
    
    // Constructor
    SynergisticInteraction(const std::string& enhancer, 
                          const std::string& enhanced, 
                          const Parameters& params);
    
    // Calculate how much target protein is enhanced by the source
    std::vector<std::pair<std::string, double>> calculateEffect(
        double sourceAmount, 
        double targetAmount, 
        double dt,
        double availableATP) const override;
    
private:
    double m_saturationConstant;  // For Hill-type kinetics
    double m_enhancementDuration; // How long the enhancement lasts
}; 