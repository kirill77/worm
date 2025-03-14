#pragma once

#include "ProteinInteraction.h"

/**
 * Represents an antagonistic interaction between proteins
 * where one protein negatively affects another (e.g., phosphorylation).
 */
class AntagonisticInteraction : public ProteinInteraction
{
public:
    // Parameters specific to antagonistic interactions
    struct Parameters {
        double removalRate;             // Rate of target protein removal
        double recoveryRate;            // Rate of target protein recovery
        double saturationConstant;      // Saturation constant for Hill-type kinetics
        ProteinInteraction::Mechanism mechanism;  // Mechanism of antagonism
        double atpCost;                 // ATP consumed per unit of protein removed
    };
    
    // Constructor
    AntagonisticInteraction(const std::string& antagonist, 
                           const std::string& target, 
                           const Parameters& params);
    
    // Calculate how much target protein is removed by the antagonist
    std::vector<std::pair<std::string, double>> calculateEffect(
        double sourceAmount, 
        double targetAmount, 
        double dt,
        double availableATP) const override;
    
    // Calculate recovery of affected protein
    double calculateRecovery(double removedAmount, double dt) const;
    
    // Get the saturation constant
    double getSaturationConstant() const { return m_saturationConstant; }
    
private:
    double m_saturationConstant;  // For Hill-type kinetics
}; 