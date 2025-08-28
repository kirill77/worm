#pragma once

#include <string>
#include <memory>
#include <vector>
#include "Molecule.h"

class MRNA
{
private:
    std::string m_sName;         // Name of gene/protein (same for mRNA)
    double m_fNumber;            // Amount of this mRNA in the cell
    double m_fHalfLife;          // How quickly it degrades (in time units)
    double m_fTranslationRate;   // Rate of protein production

public:
    MRNA(const std::string& name, double number, double halfLife, double translationRate = 1.0)
        : m_sName(name),
          m_fNumber(number), m_fHalfLife(halfLife),
          m_fTranslationRate(translationRate)
    {
        // Name represents both the gene name and the protein it produces
    }

    // Getters
    std::string getName() const { return m_sName; }
    std::string getGeneName() const { return m_sName; }    // For backward compatibility
    std::string getProteinName() const { return m_sName; } // For backward compatibility
    double getNumber() const { return m_fNumber; }
    
    // Setters
    void setNumber(double number) { m_fNumber = number; }
    void addNumber(double amount) { m_fNumber += amount; }
    
    // RNA degradation
    void degrade(double dt)
    {
        // Simple exponential decay model for mRNA degradation
        m_fNumber *= exp(-dt / m_fHalfLife);
    }

    // Translation
    std::shared_ptr<MPopulation> translate(double dt, const std::vector<std::shared_ptr<class TRNA>>& availableTRNAs) const;

    // RNA processing
    void splice();  // For future implementation of RNA processing
};

