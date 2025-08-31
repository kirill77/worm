#pragma once

#include <string>
#include <memory>
#include <assert.h>
#include <functional>
#include <cstdint>

// Forward declaration
class BindingSurface;

// Chemical type classification for molecules
enum class ChemicalType : uint8_t {
    PROTEIN,        // Any amino acid chain
    AMINO_ACID,     // Single amino acid
    DNA,           // DNA polymers
    RNA,           // RNA polymers  
    NUCLEOTIDE,    // Single nucleotides (ATP, GTP, dATP, etc.)
    LIPID,         // Fatty acids, phospholipids, steroids
    ION,           // Charged atoms/molecules (Na+, Cl-, etc.)
    OTHER          // Catch-all for everything else
};

// Simple molecule class with just name
class Molecule
{
public:
    // Default constructor
    Molecule() : m_sName(), m_type(ChemicalType::OTHER) {}
    
    // Constructor with name
    Molecule(const std::string& name) : m_sName(name), m_type(ChemicalType::OTHER) {}
    
    // Constructor with name and type
    Molecule(const std::string& name, ChemicalType type) : m_sName(name), m_type(type) {}
    
    // Get name
    const std::string& getName() const { return m_sName; }
    
    // Get chemical type
    ChemicalType getType() const { return m_type; }
    
    // Equality operator for unordered_map (based on name only for compatibility)
    bool operator==(const Molecule& other) const {
        return m_sName == other.m_sName;
    }
    
    bool operator!=(const Molecule& other) const {
        return !(*this == other);
    }

private:
    std::string m_sName;  // Name/type of the molecule
    ChemicalType m_type;  // Chemical classification of the molecule
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

// Molecule population - simple composition of Molecule + Population
struct MPopulation
{
public:
    Molecule m_molecule;     // What molecule this is
    Population m_population; // Population properties (count, binding, etc.)
    
    // Constructors
    MPopulation(const std::string& sName, double fNumber)
        : m_molecule(sName), m_population(fNumber) {}
        
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
            return hash<string>()(molecule.getName());
        }
    };
}