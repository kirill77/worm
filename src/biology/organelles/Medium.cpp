#include "pch.h"
#include "Medium.h"
#include "chemistry/MoleculeWiki.h"
#include "chemistry/ResourceDistributor.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "chemistry/GridCell.h"

// Global random number generator for consistent randomness
static std::mt19937 g_rng(std::random_device{}());

void Medium::addProtein(const MPopulation& protein, const float3& position)
{
    GridCell& gridCell = m_grid.findCell(position);
    MPopulation& cellProtein = gridCell.getOrCreateMolecule(protein.getName());

    cellProtein.bindTo(protein.getBindingSurface());
    cellProtein.m_population.m_fNumber += protein.m_population.m_fNumber;
}

void Medium::addMRNA(std::shared_ptr<MRNA> pMRNA, const float3& position)
{
    MRNA& existingMRNA = m_grid.findCell(position).getOrCreateMRNA(pMRNA->getName());
    
    // If this is a newly created mRNA (with default values), copy the properties
    if (existingMRNA.getNumber() == 0.0) {
        existingMRNA = *pMRNA;
    } else {
        // Add to existing mRNA count (accumulate molecules of same type)
        existingMRNA.addNumber(pMRNA->getNumber());
    }
}

void Medium::addTRNA(std::shared_ptr<TRNA> pTRNA, const float3& position)
{
    m_grid.findCell(position).m_pTRNAs.push_back(pTRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itProtein = gridCell.m_molecules.find(proteinName);
    return (itProtein != gridCell.m_molecules.end()) ? itProtein->second.m_population.m_fNumber : 0.0;
}

void Medium::updateProteinInteraction(double fDt)
{
    // Get all protein interactions 
    const auto& vecInteractions = MoleculeWiki::GetProteinInteractions();
    
    // First, apply direct protein interactions
    for (size_t uCell = 0; uCell < m_grid.size(); ++uCell)
    {
        // at first make a dry run to figure out who needs which resources
        m_resDistributor.notifyNewDryRun(m_grid[uCell]);
        for (size_t i = 0; i < vecInteractions.size(); ++i)
        {
            m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]);
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // now do the real run to distribute the resources
        m_resDistributor.notifyNewRealRun();
        for (size_t i = 0; i < vecInteractions.size(); ++i)
        {
            if (!m_resDistributor.notifyNewInteractionStarting(*vecInteractions[i]))
                continue;
            vecInteractions[i]->apply(m_grid[uCell], fDt, m_resDistributor);
        }

        // Ensure ATP doesn't go below zero
        auto& atpMolecule = m_grid[uCell].getOrCreateMolecule("ATP");
        atpMolecule.m_population.m_fNumber = std::max(0.0, atpMolecule.m_population.m_fNumber);
    }
}

double Medium::getTotalProteinNumber(const std::string& proteinName) const
{
    double fTotal = 0.0;
    for (uint32_t uCell = 0; uCell < m_grid.size(); ++uCell)
    {
        const GridCell& gridCell = m_grid[uCell];
        auto itProtein = gridCell.m_molecules.find(proteinName);
        if (itProtein != gridCell.m_molecules.end()) {
            fTotal += itProtein->second.m_population.m_fNumber;
        }
    }
    return fTotal;
}

void Medium::update(double fDt)
{
    m_diffusion.updateDiffusion(m_grid, fDt);
    
    // Update tRNA charging in all grid cells
    for (size_t i = 0; i < m_grid.size(); ++i) {
        m_grid[i].updateTRNAs(fDt);
    }
    
    // Interaction of proteins between each other
    updateProteinInteraction(fDt);
    
    // Update mRNA positions
    translateMRNAs(fDt);
}

void Medium::translateMRNAs(double fDt)
{
    static constexpr double ATP_PER_TRANSLATION = 4.0;  // ATP cost per amino acid (approximate)
    
    // Process translation in each grid cell
    for (size_t i = 0; i < m_grid.size(); ++i) {
        GridCell& cell = m_grid[i];
        
        // Skip cells with no mRNAs
        if (!cell.hasMRNAs()) continue;
        
        // Collect available tRNAs from this cell
        std::vector<std::shared_ptr<TRNA>> availableTRNAs;
        for (const auto& pTRNA : cell.m_pTRNAs) {
            if (pTRNA->isCharged() && pTRNA->getNumber() > 0.1) {
                availableTRNAs.push_back(pTRNA);
            }
        }
        
        // Skip if no charged tRNAs available
        if (availableTRNAs.empty()) continue;
        
        // Attempt translation for each mRNA
        auto& mrnas = cell.getMRNAs();
        auto it = mrnas.begin();
        while (it != mrnas.end()) {
            auto& mrna = it->second;
            
            // Check if we have enough ATP for translation
            // (approximate cost - in reality would depend on protein length)
            float3 cellPosition = m_grid.indexToPosition(i);
            if (getAvailableATP(cellPosition) < ATP_PER_TRANSLATION * 10) {
                ++it;
                continue;
            }
            
            // Attempt translation
            auto pProtein = mrna.translate(fDt, availableTRNAs);
            
            if (pProtein && pProtein->m_population.m_fNumber > 0.0) {
                // Translation successful - add protein to cell
                MPopulation& cellProtein = cell.getOrCreateMolecule(pProtein->getName());
                cellProtein.m_population.m_fNumber += pProtein->m_population.m_fNumber;
                
                // Consume ATP (simplified - should be proportional to protein length)
                double atpCost = ATP_PER_TRANSLATION * pProtein->m_population.m_fNumber;
                consumeATP(atpCost, cellPosition);
                
                // Reduce mRNA amount (simplified degradation from translation)
                // In reality, ribosomes can translate the same mRNA multiple times
            }
            
            ++it;
        }
    }
}

void Medium::addATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpMolecule = gridCell.getOrCreateMolecule("ATP");
    atpMolecule.m_population.m_fNumber = std::min<double>(atpMolecule.m_population.m_fNumber + fAmount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpMolecule = gridCell.getOrCreateMolecule("ATP");
    if (atpMolecule.m_population.m_fNumber >= fAmount)
    {
        atpMolecule.m_population.m_fNumber -= fAmount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itMolecule = gridCell.m_molecules.find("ATP");
    return (itMolecule != gridCell.m_molecules.end()) ? itMolecule->second.m_population.m_fNumber : 0.0;
}

Medium::Medium()
{
    m_fVolumeMicroM = 23561.0;  // Initialize volume to 23561 micrometers
    // No need to initialize protein antagonisms here anymore
    // They are now managed by MoleculeWiki
}
