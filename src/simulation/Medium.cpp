#include "pch.h"
#include "Medium.h"
#include "ProteinWiki.h"
#include "ResourceDistributor.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "GridCell.h"

// Global random number generator for consistent randomness
static std::mt19937 g_rng(std::random_device{}());

void Medium::addProtein(const MPopulation& protein, const float3& position)
{
    GridCell& gridCell = m_grid.findCell(position);
    MPopulation& cellProtein = gridCell.getOrCreateMolecule(protein.m_sName);

    cellProtein.bindTo(protein.getBindingSurface());
    cellProtein.m_fNumber += protein.m_fNumber;
}

void Medium::addMRNA(std::shared_ptr<MRNA> pMRNA, const float3& position)
{
    m_grid.findCell(position).m_pMRNAs.push_back(pMRNA);
}

double Medium::getProteinNumber(const std::string& proteinName, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itProtein = gridCell.m_molecules.find(proteinName);
    return (itProtein != gridCell.m_molecules.end()) ? itProtein->second.m_fNumber : 0.0;
}

void Medium::updateProteinInteraction(double fDt)
{
    // Get all protein interactions 
    const auto& vecInteractions = ProteinWiki::GetProteinInteractions();
    
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
        atpMolecule.m_fNumber = std::max(0.0, atpMolecule.m_fNumber);
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
            fTotal += itProtein->second.m_fNumber;
        }
    }
    return fTotal;
}

void Medium::update(double fDt)
{
    m_diffusion.updateDiffusion(m_grid, fDt);
    
    // Interaction of proteins between each other
    updateProteinInteraction(fDt);
    
    // Update mRNA positions
    translateMRNAs(fDt);
}

void Medium::translateMRNAs(double fDt)
{
    // TODO: Implement translation of mRNAs into proteins
    // This will need to:
    // 1. Check for available tRNAs
    // 2. Create new proteins
    // 3. Add proteins to appropriate cytoplasmic regions
}

void Medium::addATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpMolecule = gridCell.getOrCreateMolecule("ATP");
    atpMolecule.m_fNumber = std::min<double>(atpMolecule.m_fNumber + fAmount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpMolecule = gridCell.getOrCreateMolecule("ATP");
    if (atpMolecule.m_fNumber >= fAmount)
    {
        atpMolecule.m_fNumber -= fAmount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itMolecule = gridCell.m_molecules.find("ATP");
    return (itMolecule != gridCell.m_molecules.end()) ? itMolecule->second.m_fNumber : 0.0;
}

Medium::Medium()
{
    // No need to initialize protein antagonisms here anymore
    // They are now managed by ProteinWiki
}
