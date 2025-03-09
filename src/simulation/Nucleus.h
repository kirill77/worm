#pragma once

#include "Organelle.h"
#include "DNA.h"
#include "Chromosome.h"
#include <memory>
#include <vector>

// Forward declarations
class MRNA;

class Nucleus : public Organelle
{
private:
    std::vector<Chromosome> m_chromosomes;
    double m_fEnvelopeIntegrity;  // 1.0 = intact, 0.0 = broken down
    
    static constexpr double fENVELOPE_BREAKDOWN_RATE = 0.2f;  // Rate of nuclear envelope breakdown
    static constexpr double fENVELOPE_REFORM_RATE = 0.5f;    // Rate of nuclear envelope reformation

public:
    Nucleus(std::shared_ptr<DNA> pDNA, size_t numChromosomes = 6)  // Default to C. elegans (6 chromosomes)
        : m_fEnvelopeIntegrity(1.0f)
    {
        // Initialize chromosomes, distributing DNA among them
        m_chromosomes.reserve(numChromosomes);
        for (size_t i = 0; i < numChromosomes; ++i)
        {
            // For now, each chromosome gets a copy of the full DNA
            // In a more sophisticated simulation, we would partition the DNA
            m_chromosomes.emplace_back(pDNA);
        }
    }

    void update(double fDt, Cell& cell, Medium& medium) override;

    // Chromosome-related functions
    bool areChromosomesCondensed() const;
    bool areChromosomesAttached() const;
    bool areChromosomesSeparated() const;
    bool areChromosomesDecondensed() const;
    std::vector<std::shared_ptr<MRNA>> transcribeAll(double fDt) const;

    // Getters
    double getEnvelopeIntegrity() const { return m_fEnvelopeIntegrity; }
    size_t getChromosomeCount() const { return m_chromosomes.size(); }
    const std::vector<Chromosome>& getChromosomes() const { return m_chromosomes; }
};

