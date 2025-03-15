#include "pch.h"
#include "Membrane.h"
#include "Medium.h"
#include "log/ILog.h"
#include <algorithm>
#include <cmath>

Membrane::Membrane(std::shared_ptr<Medium> pInternalMedium, double fThickness, double fSurfaceArea)
    : ProteinBindingSurface(fSurfaceArea)  // Call base class constructor with surface area and binding capacity
    , m_pInternalMedium(pInternalMedium)
    , m_fThickness(fThickness)
{
    // Ensure the internal medium is valid
    if (!m_pInternalMedium) {
        LOG_ERROR("Internal medium cannot be null");
        // No exception throwing, but we should still return to prevent operating on null
        return;
    }
}

void Membrane::update(double dt)
{
    // Update internal medium - its dynamics are independent of external medium
    m_pInternalMedium->update(dt);
    
    // Note: In a more advanced implementation, this method could include:
    // - Membrane fluidity changes
    // - Lipid raft formation/movement
    // - Membrane protein reorganization
    // - Passive transport based on concentration gradients
    // - Signal transduction
}

bool Membrane::transportProteinInward(Medium& externalMedium, 
                                     const std::string& proteinName, 
                                     double amount, 
                                     const float3& position)
{
    // Check if the external medium has enough of the protein
    if (externalMedium.getProteinNumber(proteinName, position) < amount) {
        return false; // Not enough protein available
    }
    
    // Create a protein population for the amount we want to transport
    ProteinPopulation protein(proteinName, amount);
    
    // Remove from external medium
    // Note: This is a simplification. In reality, we would need to modify
    // the protein count directly in the external medium.
    float3 externalPos = position;
    ProteinPopulation externalProtein(proteinName, -amount); // Negative amount for removal
    externalMedium.addProtein(externalProtein, externalPos);
    
    // Add to internal medium
    float3 internalPos = position; // Same position in internal medium
    m_pInternalMedium->addProtein(protein, internalPos);
    
    return true;
}

bool Membrane::transportProteinOutward(Medium& externalMedium,
                                      const std::string& proteinName,
                                      double amount,
                                      const float3& position)
{
    // Check if the internal medium has enough of the protein
    if (m_pInternalMedium->getProteinNumber(proteinName, position) < amount) {
        return false; // Not enough protein available
    }
    
    // Create a protein population for the amount we want to transport
    ProteinPopulation protein(proteinName, amount);
    
    // Remove from internal medium
    float3 internalPos = position;
    ProteinPopulation internalProtein(proteinName, -amount); // Negative amount for removal
    m_pInternalMedium->addProtein(internalProtein, internalPos);
    
    // Add to external medium
    float3 externalPos = position; // Same position in external medium
    externalMedium.addProtein(protein, externalPos);
    
    return true;
}

bool Membrane::transportATPInward(Medium& externalMedium,
                                 double amount,
                                 const float3& position)
{
    // Check if the external medium has enough ATP
    if (externalMedium.getAvailableATP(position) < amount) {
        return false; // Not enough ATP available
    }
    
    // Consume ATP from external medium
    if (!externalMedium.consumeATP(amount, position)) {
        return false; // Could not consume ATP
    }
    
    // Add ATP to internal medium
    m_pInternalMedium->addATP(amount, position);
    
    return true;
}

bool Membrane::transportATPOutward(Medium& externalMedium,
                                  double amount,
                                  const float3& position)
{
    // Check if the internal medium has enough ATP
    if (m_pInternalMedium->getAvailableATP(position) < amount) {
        return false; // Not enough ATP available
    }
    
    // Consume ATP from internal medium
    if (!m_pInternalMedium->consumeATP(amount, position)) {
        return false; // Could not consume ATP
    }
    
    // Add ATP to external medium
    externalMedium.addATP(amount, position);
    
    return true;
} 