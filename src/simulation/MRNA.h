#pragma once

#include <string>
#include <memory>
#include <vector>
#include "Protein.h"

class MRNA
{
private:
    std::string m_sGeneName;     // Gene that produced this mRNA
    std::string m_sProteinName;  // Protein it produces
    double m_fNumber;            // Amount of this mRNA in the cell
    double m_fHalfLife;          // How quickly it degrades (in time units)
    double m_fTranslationRate;   // Rate of protein production

public:
    MRNA(const std::string& geneName, double number, double halfLife, double translationRate = 1.0)
        : m_sGeneName(geneName),
          m_fNumber(number), m_fHalfLife(halfLife),
          m_fTranslationRate(translationRate)
    {
        // Derive protein name from gene name (in practice, this would be more complex)
        m_sProteinName = geneName + "_protein";
    }

    // Getters
    std::string getGeneName() const { return m_sGeneName; }
    std::string getProteinName() const { return m_sProteinName; }
    double getNumber() const { return m_fNumber; }
    
    // RNA degradation
    void degrade(double dt)
    {
        // Simple exponential decay model for mRNA degradation
        m_fNumber *= exp(-dt / m_fHalfLife);
    }

    // Translation
    std::shared_ptr<Protein> translate(double dt, const std::vector<std::shared_ptr<class TRNA>>& availableTRNAs) const;

    // RNA processing
    void splice();  // For future implementation of RNA processing
};

