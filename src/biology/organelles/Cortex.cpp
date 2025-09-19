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
    
    // Note: TensionSphere physics is now handled by CellSim
    // Note: In a more advanced implementation, this method could include:
    // - Membrane fluidity changes
    // - Lipid raft formation/movement
    // - Membrane protein reorganization
    // - Passive transport based on concentration gradients
    // - Signal transduction
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
                
                // Create species-specific binding site protein and add to the medium
                Species species = pCell->getSpecies();
                MPopulation bindingSites(Molecule(StringDict::ID::ORGANELLE_CORTEX, ChemicalType::PROTEIN, species), amountPerPosition);
                bindingSites.bindTo(shared_from_this());
                internalMedium.addMolecule(bindingSites, normalizedPos);
            }
        }
    }
    
    return true;
} 