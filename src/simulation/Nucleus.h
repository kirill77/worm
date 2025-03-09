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
    std::shared_ptr<DNA> m_pDNA;
    std::vector<Chromosome> m_chromosomes;
    double m_envelopeIntegrity;  // 1.0 = intact, 0.0 = broken down
    
    static constexpr double ENVELOPE_BREAKDOWN_RATE = 0.2;  // Rate of nuclear envelope breakdown
    static constexpr double ENVELOPE_REFORM_RATE = 0.5;    // Rate of nuclear envelope reformation

public:
    Nucleus(std::shared_ptr<DNA> pDNA, size_t numChromosomes = 6)  // Default to C. elegans (6 chromosomes)
        : m_pDNA(pDNA)
        , m_envelopeIntegrity(1.0)
    {
        // Initialize chromosomes
        m_chromosomes.resize(numChromosomes);
    }

    void update(double dt, Cell& cell, Medium& medium) override;

    // Chromosome-related functions
    bool areChromosomesCondensed() const;
    bool areChromosomesAttached() const;
    bool areChromosomesSeparated() const;
    bool areChromosomesDecondensed() const;

    // Getters
    double getEnvelopeIntegrity() const { return m_envelopeIntegrity; }
    size_t getChromosomeCount() const { return m_chromosomes.size(); }
};

