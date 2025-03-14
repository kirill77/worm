#include "pch.h"
#include "SynergisticInteraction.h"
#include <algorithm>
#include <cmath>

SynergisticInteraction::SynergisticInteraction(
    const std::string& enhancer, 
    const std::string& enhanced, 
    const Parameters& params)
    : ProteinInteraction(enhancer, enhanced, 
                       InteractionType::SYNERGISTIC, 
                       params.mechanism, 
                       params.enhancementFactor, 
                       0.0,  // No recovery rate for synergistic interactions
                       params.atpCost)
    , m_saturationConstant(params.saturationConstant)
    , m_enhancementDuration(params.enhancementDuration)
{
}

std::vector<std::pair<std::string, double>> SynergisticInteraction::calculateEffect(
    double sourceAmount, 
    double targetAmount, 
    double dt,
    double availableATP) const
{
    std::vector<std::pair<std::string, double>> effects;
    
    // No effect if either protein is missing
    if (sourceAmount <= 0 || targetAmount <= 0) {
        return effects;
    }
    
    // Calculate enhancement using Hill-like kinetics
    double enhancementRate = m_strength * sourceAmount / (m_saturationConstant + sourceAmount);
    
    // Calculate amount enhanced in this time step
    double enhancedAmount = enhancementRate * targetAmount * dt;
    
    // Check ATP availability for energy-requiring mechanisms
    if (m_mechanism == Mechanism::BINDING || m_mechanism == Mechanism::ACTIVATION) {
        double requiredATP = enhancedAmount * m_atpCost;
        if (availableATP < requiredATP) {
            // Scale down effect if not enough ATP
            enhancedAmount = enhancedAmount * (availableATP / requiredATP);
        }
    }
    
    // Apply the effect (positive because this is enhancement)
    effects.push_back(std::make_pair(m_targetProtein, enhancedAmount));
    
    // Some mechanisms might have effects on the source protein too
    if (m_mechanism == Mechanism::BINDING) {
        // In binding, the source protein also gets "used up" in the complex
        double sourceUsed = enhancedAmount * 0.1;  // Assume 10% of source is used
        effects.push_back(std::make_pair(m_sourceProtein, -sourceUsed));
    }
    
    return effects;
} 