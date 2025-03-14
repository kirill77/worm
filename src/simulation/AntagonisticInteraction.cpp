#include "pch.h"
#include "AntagonisticInteraction.h"
#include <algorithm>
#include <cmath>

AntagonisticInteraction::AntagonisticInteraction(
    const std::string& antagonist, 
    const std::string& target, 
    const Parameters& params)
    : ProteinInteraction(antagonist, target, 
                       InteractionType::ANTAGONISTIC, 
                       params.mechanism, 
                       params.removalRate, 
                       params.recoveryRate, 
                       params.atpCost)
    , m_saturationConstant(params.saturationConstant)
{
}

std::vector<std::pair<std::string, double>> AntagonisticInteraction::calculateEffect(
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
    
    // Calculate removal using Hill-like kinetics
    double removalRate = m_strength * sourceAmount / (m_saturationConstant + sourceAmount);
    
    // Calculate amount removed in this time step
    double removedAmount = removalRate * targetAmount * dt;
    
    // Check ATP availability for phosphorylation mechanisms
    if (m_mechanism == Mechanism::PHOSPHORYLATION) {
        double requiredATP = removedAmount * m_atpCost;
        if (availableATP < requiredATP) {
            // Scale down effect if not enough ATP
            removedAmount = removedAmount * (availableATP / requiredATP);
        }
    }
    
    // Apply the effect (negative because this is removal)
    effects.push_back(std::make_pair(m_targetProtein, -removedAmount));
    
    return effects;
}

double AntagonisticInteraction::calculateRecovery(double removedAmount, double dt) const
{
    // Calculate recovery based on previously removed amount
    return removedAmount * m_recoveryRate * dt;
} 