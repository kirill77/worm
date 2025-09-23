#include "ResourceDistributor.h"
#include "GridCell.h"
#include "MoleculeInteraction.h"
#include <cassert>
#include <algorithm>

ResourceDistributor::ResourceDistributor()
{
    // Initialize with default values
}

ResourceDistributor::~ResourceDistributor()
{
    // Cleanup if needed
}

void ResourceDistributor::notifyNewDryRun(const class GridCell& cell)
{
    // Start a new dry run - increment the ID to mark a new set of computations
    ++m_curDryRunId;

    updateAvailableResources(cell);
}

bool ResourceDistributor::notifyNewInteractionStarting(const MoleculeInteraction& interaction)
{
    // Look up or create an entry for this interaction
    auto& interactionData = m_interactions[&interaction];
    m_pCurInteraction = &interactionData;
    if (isDryRun())
    {
        m_pCurInteraction->m_fScalingFactor = 1;
        m_pCurInteraction->m_requestedMolecules.resize(0);
        return true;
    }
    // previously we had real run - it had to have scaling factor of 1 or 0
    assert(m_pCurInteraction->m_fScalingFactor == 1 || m_pCurInteraction->m_fScalingFactor == 0);
    if (m_pCurInteraction->m_lastValidDryRunId != m_curDryRunId ||
        m_pCurInteraction->m_fScalingFactor == 0)
    {
        // the interaction didn't request any resources - we can skip it
        return false;
    }
    // update the scaling factor
    for (const auto &molecule : m_pCurInteraction->m_requestedMolecules)
    {
        auto it = m_resources.find(molecule);
        // if the interaction needs a resource that's not available - it can't run
        if (it == m_resources.end() || it->second.m_dryRunId != m_curDryRunId)
        {
            return false;
        }
        double fResourceScalingFactor = it->second.computeScalingFactor();
        // the interaction is constrained by the most scarce resource
        m_pCurInteraction->m_fScalingFactor = std::min(m_pCurInteraction->m_fScalingFactor,
            fResourceScalingFactor);
    }
    return true;
}

double ResourceDistributor::getAvailableResource(const Molecule& molecule)
{
    auto it = m_resources.find(molecule);
    // if we don't have such resource, or the data is stale
    if (it == m_resources.end() || it->second.m_dryRunId != m_curDryRunId)
        return 0;
    return it->second.m_fAvailable * m_pCurInteraction->m_fScalingFactor;
}

void ResourceDistributor::notifyResourceWanted(const Molecule& molecule, double amount)
{
    assert(amount > 0); // seems sub-optimal - this interaction must have bailed out earlier

    auto it = m_resources.find(molecule);
    
    // if we don't have such resource - can't distribute it
    if (it == m_resources.end())
    {
        // seems sub-optimal - this interaction must have bailed out earlier
        // if it's ATP - it's fine because if it's unavailable - it's some rare corner case
        assert(molecule == Molecule(StringDict::ID::ATP, ChemicalType::NUCLEOTIDE));
        // mark this interaction as invalid
        m_pCurInteraction->m_fScalingFactor = 0;
        return;
    }

    it->second.m_fRequested += amount;
    m_pCurInteraction->m_requestedMolecules.push_back(molecule);
    m_pCurInteraction->m_lastValidDryRunId = m_curDryRunId;
}

void ResourceDistributor::notifyNewRealRun()
{
    assert(m_curRealRunId < m_curDryRunId);
    m_curRealRunId = m_curDryRunId;
}

void ResourceDistributor::updateAvailableResources(const GridCell &cell)
{
    // Then update the available amounts
    for (const auto& [molecule, population] : cell.m_molecules)
    {
        // Update the available amount for this molecule
        auto& resource = m_resources[molecule];
        resource.m_fAvailable = population.m_fNumber;
        resource.m_dryRunId = m_curDryRunId;
    }
}
