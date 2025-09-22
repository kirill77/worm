#include "pch.h"
#include "Medium.h"
#include "chemistry/interactions/InteractionsWiki.h"
#include "chemistry/interactions/ResourceDistributor.h"
#include "chemistry/molecules/TRNA.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "chemistry/interactions/GridCell.h"

// Global random number generator for consistent randomness
static std::mt19937 g_rng(std::random_device{}());

void Medium::addMolecule(const MPopulation& population, const float3& position)
{
    GridCell& gridCell = m_grid.findCell(position);
    Population& moleculePop = gridCell.getOrCreateMolPop(population.m_molecule);

    // it's the same molecule - so they're either both bound, or both unbound
    assert(moleculePop.m_fNumber == 0.0 || moleculePop.isBound() == population.isBound());

    moleculePop.setBound(population.isBound());
    moleculePop.m_fNumber += population.m_population.m_fNumber;
}

void Medium::toBindingSites(std::vector<BindingSite>& bindingSites, const std::vector<Molecule>& bindableMolecules)
{
    // Group binding sites by grid cell index
    std::unordered_map<uint32_t, std::vector<size_t>> cellToSites;
    cellToSites.reserve(bindingSites.size());
    for (size_t i = 0; i < bindingSites.size(); ++i)
    {
        const float3& pos = bindingSites[i].m_normalized;
        uint32_t cellIndex = m_grid.positionToIndex(pos);
        cellToSites[cellIndex].push_back(i);
    }

    // For each cell, distribute bindable molecules uniformly across all binding sites in that cell
    for (const auto& entry : cellToSites)
    {
        uint32_t cellIndex = entry.first;
        const std::vector<size_t>& siteIndices = entry.second;
        if (siteIndices.empty()) continue;

        GridCell& gridCell = m_grid[cellIndex];
        const size_t numSites = siteIndices.size();

        for (const Molecule& mol : bindableMolecules)
        {
            auto it = gridCell.m_molecules.find(mol);
            if (it == gridCell.m_molecules.end())
                continue;

            Population& cellPop = it->second;
            if (cellPop.m_fNumber <= 0.0)
                continue;
            assert(cellPop.isBound()); // we're working with binding sites here

            double totalAmount = cellPop.m_fNumber;
            double share = totalAmount / static_cast<double>(numSites);

            for (size_t idx : siteIndices)
            {
                BindingSite& site = bindingSites[idx];
                auto& pop = site.m_bsMolecules[mol];
                pop.m_fNumber += share;
                pop.setBound(true); // we're working with binding sites here
            }

            // Remove from grid cell after distribution
            cellPop.m_fNumber = 0.0;
            gridCell.m_molecules.erase(it);
        }
    }
}


double Medium::getMoleculeNumber(const Molecule& molecule, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itMolecule = gridCell.m_molecules.find(molecule);
    return (itMolecule != gridCell.m_molecules.end()) ? itMolecule->second.m_fNumber : 0.0;
}

void Medium::updateMoleculeInteraction(double fDt)
{
    // Get all protein interactions 
    const auto& vecInteractions = InteractionsWiki::GetMoleculeInteractions();
    
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
        auto& atpPop = m_grid[uCell].getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
        atpPop.m_fNumber = std::max(0.0, atpPop.m_fNumber);
    }
}

double Medium::getTotalMoleculeNumber(const Molecule& molecule) const
{
    double fTotal = 0.0;
    for (uint32_t uCell = 0; uCell < m_grid.size(); ++uCell)
    {
        const GridCell& gridCell = m_grid[uCell];
        auto itMolecule = gridCell.m_molecules.find(molecule);
        if (itMolecule != gridCell.m_molecules.end()) {
            fTotal += itMolecule->second.m_fNumber;
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
    updateMoleculeInteraction(fDt);
    
    // Translation is now handled by MoleculeInteraction system
}


void Medium::addATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpPop = gridCell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
    atpPop.m_fNumber = std::min<double>(atpPop.m_fNumber + fAmount, MAX_ATP_PER_CELL);
}

bool Medium::consumeATP(double fAmount, const float3& position)
{
    auto& gridCell = m_grid.findCell(position);
    auto& atpPop = gridCell.getOrCreateMolPop(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
    if (atpPop.m_fNumber >= fAmount)
    {
        atpPop.m_fNumber -= fAmount;
        return true;
    }
    return false;
}

double Medium::getAvailableATP(const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itMolecule = gridCell.m_molecules.find(Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
    return (itMolecule != gridCell.m_molecules.end()) ? itMolecule->second.m_fNumber : 0.0;
}

Medium::Medium()
{
    m_fVolumeMicroM = 23561.0;  // Initialize volume to 23561 micrometers
    // No need to initialize protein antagonisms here anymore
    // They are now managed by MoleculeWiki
}
