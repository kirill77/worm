#include "pch.h"
#include "Medium.h"
#include "chemistry/MoleculeWiki.h"
#include "chemistry/ResourceDistributor.h"
#include "chemistry/TRNA.h"
#include <random>
#include <algorithm>
#include <cassert>
#include "chemistry/GridCell.h"

// Global random number generator for consistent randomness
static std::mt19937 g_rng(std::random_device{}());

void Medium::addMolecule(const MPopulation& population, const float3& position)
{
    GridCell& gridCell = m_grid.findCell(position);
    Population& moleculePop = gridCell.getOrCreateMolPop(population.m_molecule);

    moleculePop.bindTo(population.getBindingSurface());
    moleculePop.m_fNumber += population.m_population.m_fNumber;
}


double Medium::getMoleculeNumber(const Molecule& molecule, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto itMolecule = gridCell.m_molecules.find(molecule);
    return (itMolecule != gridCell.m_molecules.end()) ? itMolecule->second.m_fNumber : 0.0;
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
    updateProteinInteraction(fDt);
    
    // Update mRNA translation
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
        
        // Attempt translation for each mRNA molecule
        auto& molecules = cell.m_molecules;
        auto it = molecules.begin();
        while (it != molecules.end()) {
            if (it->first.getType() == ChemicalType::MRNA && it->second.m_fNumber > 0.1) {
                // Check if we have enough ATP for translation
                float3 cellPosition = m_grid.indexToPosition(i);
                if (getAvailableATP(cellPosition) < ATP_PER_TRANSLATION * 10) {
                    ++it;
                    continue;
                }
                
                // Get translation rate from MoleculeWiki
                const auto& info = MoleculeWiki::getInfo(it->first);
                double translationRate = info.m_fTranslationRate;
                
                // Attempt translation
                auto pProtein = it->first.translate(fDt, it->second.m_fNumber, translationRate, cell.m_molecules);
                
                if (pProtein && pProtein->m_population.m_fNumber > 0.0) {
                    // Translation successful - add protein to cell
                    Population& cellProteinPop = cell.getOrCreateMolPop(pProtein->m_molecule);
                    cellProteinPop.m_fNumber += pProtein->m_population.m_fNumber;
                    
                    // Consume ATP (simplified - should be proportional to protein length)
                    double atpCost = ATP_PER_TRANSLATION * pProtein->m_population.m_fNumber;
                    consumeATP(atpCost, cellPosition);
                    
                    // Reduce mRNA amount (simplified degradation from translation)
                    // In reality, ribosomes can translate the same mRNA multiple times
                }
            }
            ++it;
        }
    }
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
