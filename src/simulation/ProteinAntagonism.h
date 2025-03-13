#pragma once

#include <string>
#include <memory>
#include "math/vector.h"

// Class to handle antagonistic interactions between pairs of proteins
class ProteinAntagonism
{
public:
    enum class Mechanism {
        PHOSPHORYLATION,      // Kinase-mediated phosphorylation (requires ATP)
        CORTICAL_EXCLUSION,   // Competition for cortical binding sites
        RECRUITMENT,          // Recruitment of antagonistic factors
        COMPLEX_FORMATION     // Sequestration through complex formation
    };

    struct Parameters {
        double fRemovalRate;           // Rate at which antagonist affects target
        double fRecoveryRate;          // Rate at which affected proteins recover
        double fSaturationConstant;    // Constant for saturable antagonistic effect
        Mechanism eMechanism;          // Mechanism of antagonism
        double fATPCost;               // ATP cost per unit of antagonism (0 for non-ATP mechanisms)
    };

    ProteinAntagonism(
        const std::string& sAntagonist,
        const std::string& sTarget,
        const Parameters& params
    ) : m_sAntagonist(sAntagonist)
      , m_sTarget(sTarget)
      , m_params(params)
    {}

    // Calculate how much of target protein is removed based on antagonist levels
    double calculateRemoval(double fTargetAmount, double fAntagonistAmount, double fDt, double fAvailableATP = 0.0) const
    {
        double antagonistStrength = fAntagonistAmount / (fAntagonistAmount + m_params.fSaturationConstant);
        double removalAmount = fTargetAmount * antagonistStrength * m_params.fRemovalRate * fDt;
        
        // If mechanism requires ATP, check availability
        if (m_params.eMechanism == Mechanism::PHOSPHORYLATION) {
            double atpNeeded = removalAmount * m_params.fATPCost;
            if (atpNeeded > fAvailableATP) {
                removalAmount = (fAvailableATP / m_params.fATPCost);
            }
            return removalAmount;
        }
        
        return removalAmount;
    }

    // Calculate recovery of affected proteins
    double calculateRecovery(double fRemovedAmount, double fDt) const
    {
        return fRemovedAmount * m_params.fRecoveryRate * fDt;
    }

    // Getters
    const std::string& getAntagonist() const { return m_sAntagonist; }
    const std::string& getTarget() const { return m_sTarget; }
    Mechanism getMechanism() const { return m_params.eMechanism; }
    double getATPCost() const { return m_params.fATPCost; }

private:
    std::string m_sAntagonist;         // Protein doing the antagonizing
    std::string m_sTarget;             // Protein being antagonized
    Parameters m_params;               // Interaction parameters
}; 