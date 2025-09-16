#pragma once

#include <string>
#include <memory>
#include <assert.h>
#include <functional>
#include <cstdint>
#include "StringDict.h"

// Forward declarations
class BindingSurface;
struct MPopulation;

// Chemical type classification for molecules
enum class ChemicalType : uint8_t {
    PROTEIN,        // Any amino acid chain
    AMINO_ACID,     // Single amino acid
    DNA,           // DNA polymers
    MRNA,          // Messenger RNA - carries genetic information for translation
    TRNA,          // Transfer RNA - carries amino acids during translation  
    RRNA,          // Ribosomal RNA - structural component of ribosomes
    NUCLEOTIDE,    // Single nucleotides (ATP, GTP, dATP, etc.)
    LIPID,         // Fatty acids, phospholipids, steroids
    ION,           // Charged atoms/molecules (Na+, Cl-, etc.)
    OTHER          // Catch-all for everything else
};

// Biological species/organism the molecule belongs to
enum class Species : uint8_t {
    HUMAN = 0,
    C_ELEGANS = 1
};

// Population class - handles population properties without molecule identity
class Population
{
public:
    double m_fNumber;     // Number of molecules in this population
    
    // Constructors
    Population() : m_fNumber(0.0) {}
    Population(double fNumber) : m_fNumber(fNumber) {}
    
    // Check if this population is bound to any surface
    bool isBound() const
    {
        return !m_pBindingSurface.expired();
    }

    // Get a shared_ptr to the binding surface (if bound)
    std::shared_ptr<BindingSurface> getBindingSurface() const {
        return m_pBindingSurface.lock();
    }

    // Bind this population to a surface
    void bindTo(std::shared_ptr<BindingSurface> pSurface)
    {
        // if you want to bind to a different surface - have to unbind first
        assert(m_pBindingSurface.expired() || m_pBindingSurface.lock() == pSurface);
        m_pBindingSurface = pSurface;
    }

    // Unbind this population from its current surface
    void unbind()
    {
        m_pBindingSurface.reset();
    }

private:
    // Weak pointer to the surface this population is bound to (if any)
    std::weak_ptr<BindingSurface> m_pBindingSurface;
};

// Simple molecule class with optimized storage
class Molecule
{
public:
    // Default constructor
    Molecule() : m_id(StringDict::ID::eUNKNOWN), m_type(ChemicalType::OTHER), m_species(Species::HUMAN) {}
    
    // Constructor with StringDict ID and type  
    Molecule(StringDict::ID id, ChemicalType type) : m_id(id), m_type(type), m_species(Species::HUMAN) {}

    // Constructor with explicit species
    Molecule(StringDict::ID id, ChemicalType type, Species species) : m_id(id), m_type(type), m_species(species) {}
    
    // Get name from StringDict
    const std::string& getName() const {
        return StringDict::idToString(m_id);
    }
    
    // Get chemical type
    ChemicalType getType() const { return m_type; }
    
    // Get biological species
    Species getSpecies() const { return m_species; }
    
    // Get StringDict ID (for optimization purposes)
    StringDict::ID getID() const { return m_id; }
    
    // Equality operator for unordered_map
    bool operator==(const Molecule& other) const {
        return m_id == other.m_id && m_type == other.m_type;
    }
    
    bool operator!=(const Molecule& other) const {
        return !(*this == other);
    }
    

private:
    StringDict::ID m_id;   // StringDict ID  
    ChemicalType m_type;   // Chemical classification of the molecule
    Species m_species;     // Organism/species classification
};

// Molecule population - simple composition of Molecule + Population
struct MPopulation
{
public:
    Molecule m_molecule;     // What molecule this is
    Population m_population; // Population properties (count, binding, etc.)
    
    // Constructor - requires explicit Molecule object
    MPopulation(const Molecule& molecule, double fNumber)
        : m_molecule(molecule), m_population(fNumber) {}
    
    // Convenience methods that delegate to the appropriate member
    const std::string& getName() const { return m_molecule.getName(); }
    
    // Delegate population methods for convenience
    bool isBound() const { return m_population.isBound(); }
    std::shared_ptr<BindingSurface> getBindingSurface() const { return m_population.getBindingSurface(); }
    void bindTo(std::shared_ptr<BindingSurface> pSurface) { m_population.bindTo(pSurface); }
    void unbind() { m_population.unbind(); }
};

// Hash function specialization for std::unordered_map<Molecule, T>
namespace std {
    template<>
    struct hash<Molecule> {
        size_t operator()(const Molecule& molecule) const {
            // Hash based on ID and type
            size_t h1 = hash<int>()(static_cast<int>(molecule.getID()));
            size_t h2 = hash<int>()(static_cast<int>(molecule.getType()));
            return h1 ^ (h2 << 1);
        }
    };
}