#include "pch.h"
#include "Cortex.h"
#include "Medium.h"
#include "Cell.h"
#include "utils/log/ILog.h"
#include <algorithm>
#include <cmath>

Cortex::Cortex(std::weak_ptr<Cell> pCell, double fThickness)
    : Organelle(pCell)
    , m_fThickness(fThickness)
{
    // Set the binding surface type for cortex
    m_surfaceType = StringDict::ID::ORGANELLE_CORTEX;
}

void Cortex::update(double fDtSec, Cell& cell)
{
    // Update internal medium - its dynamics are independent of external medium
    cell.getInternalMedium().update(fDtSec);

    m_tensionSphere.makeTimeStep(fDtSec);
    
    // Note: In a more advanced implementation, this method could include:
    // - Membrane fluidity changes
    // - Lipid raft formation/movement
    // - Membrane protein reorganization
    // - Passive transport based on concentration gradients
    // - Signal transduction
}

bool Cortex::transportProteinInward(Medium& externalMedium, 
                                     const std::string& proteinName, 
                                     double amount, 
                                     const float3& position)
{
    // Get the cell and its internal medium
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot transport protein: cell reference is invalid");
        return false;
    }
    
    Medium& internalMedium = pCell->getInternalMedium();

    // Check if the external medium has enough of the protein
    if (externalMedium.getProteinNumber(proteinName, position) < amount) {
        return false; // Not enough protein available
    }
    
    // Create a protein population for the amount we want to transport
    MPopulation protein(proteinName, amount);
    
    // Remove from external medium
    // Note: This is a simplification. In reality, we would need to modify
    // the protein count directly in the external medium.
    float3 externalPos = position;
    MPopulation externalProtein(proteinName, -amount); // Negative amount for removal
    externalMedium.addProtein(externalProtein, externalPos);
    
    // Add to internal medium
    float3 internalPos = position; // Same position in internal medium
    internalMedium.addProtein(protein, internalPos);
    
    return true;
}

bool Cortex::transportProteinOutward(Medium& externalMedium,
                                      const std::string& proteinName,
                                      double amount,
                                      const float3& position)
{
    // Get the cell and its internal medium
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot transport protein: cell reference is invalid");
        return false;
    }
    
    Medium& internalMedium = pCell->getInternalMedium();

    // Check if the internal medium has enough of the protein
    if (internalMedium.getProteinNumber(proteinName, position) < amount) {
        return false; // Not enough protein available
    }
    
    // Create a protein population for the amount we want to transport
    MPopulation protein(proteinName, amount);
    
    // Remove from internal medium
    float3 internalPos = position;
    MPopulation internalProtein(proteinName, -amount); // Negative amount for removal
    internalMedium.addProtein(internalProtein, internalPos);
    
    // Add to external medium
    float3 externalPos = position; // Same position in external medium
    externalMedium.addProtein(protein, externalPos);
    
    return true;
}

bool Cortex::transportATPInward(Medium& externalMedium,
                                 double amount,
                                 const float3& position)
{
    // Get the cell and its internal medium
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot transport ATP: cell reference is invalid");
        return false;
    }
    
    Medium& internalMedium = pCell->getInternalMedium();

    // Check if the external medium has enough ATP
    if (externalMedium.getAvailableATP(position) < amount) {
        return false; // Not enough ATP available
    }
    
    // Consume ATP from external medium
    if (!externalMedium.consumeATP(amount, position)) {
        return false; // Could not consume ATP
    }
    
    // Add ATP to internal medium
    internalMedium.addATP(amount, position);
    
    return true;
}

bool Cortex::transportATPOutward(Medium& externalMedium,
                                  double amount,
                                  const float3& position)
{
    // Get the cell and its internal medium
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot transport ATP: cell reference is invalid");
        return false;
    }
    
    Medium& internalMedium = pCell->getInternalMedium();

    // Check if the internal medium has enough ATP
    if (internalMedium.getAvailableATP(position) < amount) {
        return false; // Not enough ATP available
    }
    
    // Consume ATP from internal medium
    if (!internalMedium.consumeATP(amount, position)) {
        return false; // Could not consume ATP
    }
    
    // Add ATP to external medium
    externalMedium.addATP(amount, position);
    
    return true;
}

bool Cortex::initializeBindingSites(double totalAmount)
{
    // Get the cell and its internal medium
    auto pCell = getCell();
    if (!pCell) {
        LOG_ERROR("Cannot initialize binding sites: cell reference is invalid");
        return false;
    }
    
    Medium& internalMedium = pCell->getInternalMedium();

    // Number of sample points to use (adjust for desired density)
    const int sampleCount = 20; 
    
    // Calculate total number of positions
    int totalPositions = sampleCount * sampleCount * sampleCount;
    
    // Calculate amount per position to distribute totalAmount evenly
    double amountPerPosition = totalAmount / totalPositions;
    
    // Add binding sites at each position in the grid
    for (int x = -sampleCount/2; x < sampleCount/2; x++) {
        for (int y = -sampleCount/2; y < sampleCount/2; y++) {
            for (int z = -sampleCount/2; z < sampleCount/2; z++) {
                // Create a normalized position between -1 and 1
                float3 normalizedPos(
                    (float)x / (sampleCount/2),
                    (float)y / (sampleCount/2),
                    (float)z / (sampleCount/2)
                );
                
                // Create binding site protein and add to the medium
                MPopulation bindingSites(StringDict::idToString(StringDict::ID::ORGANELLE_CORTEX), amountPerPosition);
                bindingSites.bindTo(shared_from_this());
                internalMedium.addProtein(bindingSites, normalizedPos);
            }
        }
    }
    
    return true;
} 