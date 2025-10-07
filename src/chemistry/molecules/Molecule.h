#pragma once

#include <string>
#include <memory>
#include <assert.h>
#include <functional>
#include <cstdint>
#include "StringDict.h"

// Forward declarations
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
    GENERIC = 0,
    C_ELEGANS = 1,
    COUNT
};

// Population class - handles population properties without molecule identity
class Population
{
public:
    double m_fNumber;     // Number of molecules in this population
    bool m_isBound = false; // Whether this population is bound to some surface
    
    // Constructors
    Population() : m_fNumber(0.0) {}
    Population(double fNumber) : m_fNumber(fNumber) {}
    
    // Check if this population is bound to any surface
    bool isBound() const { return m_isBound; }

    // Set binding state
    void setBound(bool bound) { m_isBound = bound; }

private:
};

// Simple molecule class with optimized storage
class Molecule
{
public:
    // Default constructor
    Molecule() : m_id(StringDict::ID::eUNKNOWN), m_type(ChemicalType::OTHER), m_species(Species::GENERIC) {}
    
    // Constructor with StringDict ID and type  
    Molecule(StringDict::ID id, ChemicalType type) : m_id(id), m_type(type), m_species(Species::GENERIC) {}

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
    
    // Equality operator for unordered_map (species-aware)
    bool operator==(const Molecule& other) const {
        return m_id == other.m_id && m_type == other.m_type && m_species == other.m_species;
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
    
    // Constructor - initialize from an existing Population (preserves bound state)
    MPopulation(const Molecule& molecule, const Population& population)
        : m_molecule(molecule), m_population(population) {}
    
    // Convenience methods that delegate to the appropriate member
    const std::string& getName() const { return m_molecule.getName(); }
    
    // Delegate population methods for convenience
    bool isBound() const { return m_population.isBound(); }
    void setBound(bool bound) { m_population.setBound(bound); }
};

// Hash function specialization for std::unordered_map<Molecule, T>
namespace std {
    template<>
    struct hash<Molecule> {
        size_t operator()(const Molecule& molecule) const {
            // Hash based on ID, type, and species
            size_t h1 = hash<int>()(static_cast<int>(molecule.getID()));
            size_t h2 = hash<int>()(static_cast<int>(molecule.getType()));
            size_t h3 = hash<int>()(static_cast<int>(molecule.getSpecies()));
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };
}