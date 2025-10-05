#include "pch.h"
#include "Medium.h"
#include "chemistry/interactions/InteractionsWiki.h"
#include "chemistry/interactions/ResourceDistributor.h"
#include "chemistry/molecules/TRNA.h"
#include <random>
#include <algorithm>
#include <cassert>
#include <vector>
#include "chemistry/interactions/GridCell.h"
// Use forward-declared Cortex; include header only where needed
#include "Cortex.h"

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

double Medium::getMoleculeConcentration(const Molecule& molecule, const float3& position) const
{
    const auto& gridCell = m_grid.findCell(position);
    auto it = gridCell.m_molecules.find(molecule);
    double count = (it != gridCell.m_molecules.end()) ? it->second.m_fNumber : 0.0;
    double vol = gridCell.getVolumeMicroM3();
    if (vol <= 0.0)
        return 0.0;
    return count / vol; // molecules per µm^3
}

void Medium::updateGridCellVolumes(Cortex& cortex)
{
    // Precompute world positions of all grid vertices and reuse across cells
    const uint32_t res = Grid::resolution();
    const uint32_t vres = res + 1;

    auto edgeCoord = [&](uint32_t i) { return -1.0 + 2.0 * (static_cast<double>(i) / static_cast<double>(res)); };

    // Precompute normalized coordinates for vertices along each axis
    std::vector<double> edges(vres);
    for (uint32_t i = 0; i < vres; ++i) {
        edges[i] = -1.0 + 2.0 * (static_cast<double>(i) / static_cast<double>(res));
    }

    // Precompute world positions for each vertex (ix,iy,iz) with ix,iy,iz in [0..res]
    const uint32_t vertCount = vres * vres * vres;
    std::vector<float3> worldVerts(vertCount);
    auto vindex = [&](uint32_t ix, uint32_t iy, uint32_t iz) {
        return ix * vres * vres + iy * vres + iz;
    };
    for (uint32_t ix = 0; ix < vres; ++ix)
    for (uint32_t iy = 0; iy < vres; ++iy)
    for (uint32_t iz = 0; iz < vres; ++iz)
    {
        float3 npos((float)edges[ix], (float)edges[iy], (float)edges[iz]);
        worldVerts[vindex(ix,iy,iz)] = cortex.normalizedToWorld(npos);
    }

    // Helper to compute volume of a tetrahedron
    auto tetVolume = [](const float3& a, const float3& b, const float3& c, const float3& d) {
        float3 ab = b - a, ac = c - a, ad = d - a;
        float v = dot(ab, cross(ac, ad));
        return std::abs(v) / 6.0f;
    };

    // Now iterate over cells and use precomputed vertices
    double totalGridVolume = 0.0;
    for (uint32_t ix = 0; ix < res; ++ix)
    for (uint32_t iy = 0; iy < res; ++iy)
    for (uint32_t iz = 0; iz < res; ++iz)
    {
        const float3& c000 = worldVerts[vindex(ix,   iy,   iz  )];
        const float3& c100 = worldVerts[vindex(ix+1, iy,   iz  )];
        const float3& c010 = worldVerts[vindex(ix,   iy+1, iz  )];
        const float3& c110 = worldVerts[vindex(ix+1, iy+1, iz  )];
        const float3& c001 = worldVerts[vindex(ix,   iy,   iz+1)];
        const float3& c101 = worldVerts[vindex(ix+1, iy,   iz+1)];
        const float3& c011 = worldVerts[vindex(ix,   iy+1, iz+1)];
        const float3& c111 = worldVerts[vindex(ix+1, iy+1, iz+1)];

        double vol = 0.0;
        vol += tetVolume(c000, c100, c010, c001);
        vol += tetVolume(c100, c110, c010, c111);
        vol += tetVolume(c100, c010, c001, c111);
        vol += tetVolume(c010, c001, c011, c111);
        vol += tetVolume(c100, c001, c101, c111);

        // Store volume in corresponding GridCell (center at midpoints)
        float x0 = (float)edges[ix], x1 = (float)edges[ix+1];
        float y0 = (float)edges[iy], y1 = (float)edges[iy+1];
        float z0 = (float)edges[iz], z1 = (float)edges[iz+1];
        float3 centerNorm(
            (x0 + x1) * 0.5f,
            (y0 + y1) * 0.5f,
            (z0 + z1) * 0.5f);
        GridCell& gc = m_grid.findCell(centerNorm);
        gc.setVolumeMicroM3(vol);

        totalGridVolume += vol;
    }

    // Compare total grid volume to medium volume; assert they are reasonably close
    if (m_fVolumeMicroM > 0.0)
    {
        double relError = std::abs(totalGridVolume - m_fVolumeMicroM) / m_fVolumeMicroM;
        // Allow generous tolerance due to coarse hexahedron-to-tetrahedra approximation and mapping
        assert(relError < 0.25);
    }
}

void Medium::toBindingSites(std::vector<CortexMolecules>& bindingSites, const std::vector<Molecule>& bindableMolecules)
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
                CortexMolecules& site = bindingSites[idx];
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

// Removed getTotalMoleculeNumber to encourage concentration-based logic

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
