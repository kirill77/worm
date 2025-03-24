#include "pch.h"
#include "ResourceDistributor.h"
#include "GridCell.h"
#include "ProteinInteraction.h"
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

void ResourceDistributor::notifyNewInteractionStarting(const ProteinInteraction& interaction)
{
    // Look up or create an entry for this interaction
    auto& interactionData = m_interactions[&interaction];
    m_pCurInteraction = &interactionData;
}

double ResourceDistributor::getAvailableResource(const std::string& resourceName)
{
    auto it = m_resources.find(resourceName);
    if (it == m_resources.end())
        return 0;
    return it->second.m_fAvailable;
}

void ResourceDistributor::notifyResourceWanted(const std::string& resourceName, double amount)
{
    assert(amount > 0); // seems sub-optimal - this interaction must have bailed out earlier

    auto it = m_resources.find(resourceName);
    
    // if we don't have such resource - can't distribute it
    if (it == m_resources.end())
    {
        assert(false); // seems sub-optimal - this interaction must have bailed out earlier
        return;
    }

    it->second.m_fRequested += amount;
    m_pCurInteraction->m_consumedResourceNames.push_back(resourceName);
}

void ResourceDistributor::notifyNewRealRun()
{
    assert(m_curRealRunId < m_curDryRunId);
    m_curRealRunId = m_curDryRunId;
}

void ResourceDistributor::updateAvailableResources(const GridCell &cell)
{
    // Special case for ATP - it's stored directly in GridCell
    auto& atpResource = m_resources["ATP"];
    atpResource.m_fAvailable = cell.m_fAtp;
    atpResource.m_dryRunId = m_curDryRunId;

    // Update the available amounts from the current cell state by iterating over all proteins
    for (const auto& [proteinName, proteinPop] : cell.m_proteins)
    {
        // Update the available amount for this protein
        auto& resource = m_resources[proteinName];
        resource.m_fAvailable = proteinPop.m_fNumber;
        resource.m_dryRunId = m_curDryRunId;
    }
}
