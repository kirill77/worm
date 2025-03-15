#pragma once

#include <memory>
#include <unordered_map>
#include <string>
#include "Protein.h"
#include "math/vector.h"

/**
 * Abstract base class representing any cellular structure that can bind proteins.
 * This includes membranes, organelles, and other cellular structures.
 */
class ProteinBindingSurface : public std::enable_shared_from_this<ProteinBindingSurface>
{
protected:
    // Surface properties
    double m_fSurfaceArea;      // Surface area in square micrometers
    double m_fBindingCapacity;  // Maximum protein binding capacity per square micrometer
    
    // Map of bound proteins and their quantities
    std::unordered_map<std::string, std::shared_ptr<ProteinPopulation>> m_boundProteins;

public:
    /**
     * Constructor that initializes the binding surface
     * 
     * @param fSurfaceArea Surface area in square micrometers
     * @param fBindingCapacity Maximum protein binding capacity per square micrometer
     */
    ProteinBindingSurface(double fSurfaceArea = 1.0, double fBindingCapacity = 1e6)
        : m_fSurfaceArea(fSurfaceArea)
        , m_fBindingCapacity(fBindingCapacity)
    {}
    
    /**
     * Virtual destructor for proper cleanup in derived classes
     */
    virtual ~ProteinBindingSurface() = default;
    
    /**
     * Binds a protein population to this surface
     * 
     * @param proteinName Name of the protein to bind
     * @param amount Amount of protein to bind
     * @param position Position on the surface (if applicable)
     * @return Amount of protein actually bound (may be less if at capacity)
     */
    virtual double bindProtein(const std::string& proteinName, double amount, const float3& position = float3(0,0,0))
    {
        // Create a new protein population or update existing one
        std::shared_ptr<ProteinPopulation> pProtein;
        auto it = m_boundProteins.find(proteinName);
        
        if (it != m_boundProteins.end()) {
            // Update existing population
            pProtein = it->second;
            pProtein->m_fNumber += amount;
        } else {
            // Create new population
            pProtein = std::make_shared<ProteinPopulation>(proteinName, amount);
            m_boundProteins[proteinName] = pProtein;
        }
        
        // Set the binding relationship
        pProtein->bindTo(shared_from_this());
        
        return amount;
    }
    
    /**
     * Unbinds a protein from the surface
     * 
     * @param proteinName Name of the protein to unbind
     * @param amount Amount of protein to unbind (or all if larger than bound amount)
     * @param position Position on the surface (if applicable)
     * @return Amount of protein actually unbound
     */
    virtual double unbindProtein(const std::string& proteinName, double amount, const float3& position = float3(0,0,0))
    {
        auto it = m_boundProteins.find(proteinName);
        if (it == m_boundProteins.end()) {
            return 0.0;  // No protein bound
        }
        
        std::shared_ptr<ProteinPopulation> pProtein = it->second;
        double availableAmount = pProtein->m_fNumber;
        double unboundAmount = std::min(amount, availableAmount);
        
        // Update protein population
        pProtein->m_fNumber -= unboundAmount;
        
        // Remove if completely unbound
        if (pProtein->m_fNumber <= 0) {
            pProtein->unbind();
            m_boundProteins.erase(it);
        }
        
        return unboundAmount;
    }
    
    /**
     * Gets the amount of a specific protein bound to the surface
     * 
     * @param proteinName Name of the protein to check
     * @return Amount of protein bound
     */
    virtual double getBoundProteinAmount(const std::string& proteinName) const
    {
        auto it = m_boundProteins.find(proteinName);
        return (it != m_boundProteins.end()) ? it->second->m_fNumber : 0.0;
    }
    
    /**
     * Gets the total amount of all proteins bound to the surface
     * 
     * @return Total amount of bound proteins
     */
    virtual double getTotalBoundProteinAmount() const
    {
        double total = 0.0;
        for (const auto& [_, pProtein] : m_boundProteins) {
            total += pProtein->m_fNumber;
        }
        return total;
    }
    
    /**
     * Gets the maximum binding capacity of the surface
     * 
     * @return Maximum binding capacity
     */
    virtual double getMaxCapacity() const
    {
        return m_fSurfaceArea * m_fBindingCapacity;
    }
    
    /**
     * Gets the surface area
     * 
     * @return Surface area in square micrometers
     */
    virtual double getSurfaceArea() const
    {
        return m_fSurfaceArea;
    }
    
    /**
     * Sets the surface area
     * 
     * @param fSurfaceArea New surface area in square micrometers
     */
    virtual void setSurfaceArea(double fSurfaceArea)
    {
        m_fSurfaceArea = fSurfaceArea;
    }
}; 