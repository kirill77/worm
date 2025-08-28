#pragma once

#include <memory>
#include <vector>
#include <string>
#include <unordered_map>
#include "Molecule.h"
#include "MRNA.h"
#include "TRNA.h"

// A single cell in the 3D grid representing the simulation space
class GridCell 
{
public:
    // Minimum possible resource level (to check with assertions)
    static constexpr double MIN_RESOURCE_LEVEL = 0.0;

    std::unordered_map<std::string, MPopulation> m_molecules;
    std::vector<std::shared_ptr<TRNA>> m_pTRNAs;
    
    // Constructor
    GridCell();
    
    // Helper to get or create molecule population
    MPopulation& getOrCreateMolecule(const std::string& sMoleculeName);
    
    // Helper to get or create mRNA
    MRNA& getOrCreateMRNA(const std::string& sName);
    
    // Remove mRNA by name
    void removeMRNA(const std::string& sName) { m_pMRNAs.erase(sName); }
    
    // Check if mRNAs container is empty
    bool hasMRNAs() const { return !m_pMRNAs.empty(); }
    
    // Get mRNA count
    size_t getMRNACount() const { return m_pMRNAs.size(); }
    
    // Get const reference to mRNAs for read-only access
    const std::unordered_map<std::string, MRNA>& getMRNAs() const { return m_pMRNAs; }
    
    // Get non-const reference to mRNAs for modification access
    std::unordered_map<std::string, MRNA>& getMRNAs() { return m_pMRNAs; }
    
    // mRNA management
    void updateMRNAs(double dt);  // Handle mRNA degradation and cleanup
    
    // tRNA management  
    void updateTRNAs(double dt);  // Handle tRNA charging and cleanup

private:
    std::unordered_map<std::string, MRNA> m_pMRNAs;
}; 