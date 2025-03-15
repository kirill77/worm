#pragma once

#include <string>
#include <memory>

// Forward declaration
class ProteinBindingSurface;

struct Protein
{
    std::string m_sName;  // Name/type of the protein
    
    Protein(const std::string& name) : m_sName(name) {}
};

struct ProteinPopulation : public Protein
{
    double m_fNumber; // Number of molecules in this population
    
    // Weak pointer to the surface this protein population is bound to (if any)
    std::weak_ptr<ProteinBindingSurface> m_pBindingSurface;
    
    ProteinPopulation(const std::string& sName, double fNumber)
        : Protein(sName), m_fNumber(fNumber) {}
    
    // Check if this protein population is bound to any surface
    bool isBound() const {
        return !m_pBindingSurface.expired();
    }
    
    // Get a shared_ptr to the binding surface (if bound)
    std::shared_ptr<ProteinBindingSurface> getBindingSurface() const {
        return m_pBindingSurface.lock();
    }
    
    // Bind this protein population to a surface
    void bindTo(std::shared_ptr<ProteinBindingSurface> pSurface) {
        m_pBindingSurface = pSurface;
    }
    
    // Unbind this protein population from its current surface
    void unbind() {
        m_pBindingSurface.reset();
    }
};
