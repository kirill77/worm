#pragma once

#include "Organelle.h"
#include "chemistry/DNA.h"
#include "Chromosome.h"
#include "chemistry/GridCell.h"
#include <memory>
#include <vector>

// Forward declarations
class MRNA;

class Nucleus : public Organelle
{
private:
    std::vector<Chromosome> m_chromosomes;
    double m_fEnvelopeIntegrity;  // 1.0 = intact, 0.0 = broken down
    GridCell m_nuclearCompartment;  // Nuclear chemistry compartment
    
    static constexpr double fENVELOPE_BREAKDOWN_RATE = 0.2f;  // Rate of nuclear envelope breakdown
    static constexpr double fENVELOPE_REFORM_RATE = 0.5f;    // Rate of nuclear envelope reformation

public:
    // Constructor that takes a vector of chromosomes
    Nucleus(std::weak_ptr<Cell> pCell, const std::vector<Chromosome>& chromosomes)
        : Organelle(pCell)
        , m_chromosomes(chromosomes)
        , m_fEnvelopeIntegrity(1.0f)
    {}

    void update(double fDt, Cell& cell) override;

    // Chromosome-related functions
    bool areChromosomesCondensed() const;
    bool areChromosomesAttached() const;
    bool areChromosomesSeparated() const;
    bool areChromosomesDecondensed() const;
    std::vector<std::shared_ptr<MRNA>> transcribeAll(double fDt) const;

    // Nuclear transport
    void importProtein(const std::string& proteinName, double amount);
    void exportMRNA(std::shared_ptr<MRNA> mRNA);
    
    // Nuclear compartment access
    const GridCell& getNuclearCompartment() const { return m_nuclearCompartment; }
    GridCell& getNuclearCompartment() { return m_nuclearCompartment; }
    
    // Getters
    double getEnvelopeIntegrity() const { return m_fEnvelopeIntegrity; }
    size_t getChromosomeCount() const { return m_chromosomes.size(); }
    const std::vector<Chromosome>& getChromosomes() const { return m_chromosomes; }
};

