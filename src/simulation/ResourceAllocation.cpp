#include "pch.h"
#include "ResourceAllocation.h"
#include "GridCell.h"
#include "ProteinInteraction.h"
#include <cassert>
#include <algorithm>

ResourceAllocation::ResourceAllocation()
{
    // Initialize with default values
}

ResourceAllocation::~ResourceAllocation()
{
    // Cleanup if needed
}

void ResourceAllocation::notifyNewDryRun(const class GridCell& cell)
{
    // Start a new dry run - increment the ID to mark a new set of computations
    ++m_curDryRunId;

    updateAvailableResources(cell);
}

double ResourceAllocation::notifyNewInteractionStarting(const ProteinInteraction& interaction)
{
    // Look up or create an entry for this interaction
    auto& interactionData = m_interactions[&interaction];

    // Check if this interaction has already been processed in the current dry run
    if (interactionData.m_dryRunId != m_curDryRunId)
    {
        // This is a new interaction for this dry run
        interactionData.m_dryRunId = m_curDryRunId;
        interactionData.m_consumedResourceNames.clear();
        interactionData.m_fScalingFactor = 1.0;
    }
    else
    {
        // For the real run, compute the scaling factor based on resources consumed during the dry run

        // Look through all resources this interaction consumed and find the most constraining one
        for (const auto& resourceName : interactionData.m_consumedResourceNames)
        {
            auto resourceIt = m_resources.find(resourceName);
            if (resourceIt != m_resources.end())
            {
                double resourceScaling = resourceIt->second.computeScalingFactor();

                // Take the smallest scaling factor (most constraining resource)
                interactionData.m_fScalingFactor = std::min(interactionData.m_fScalingFactor, resourceScaling);
            }
        }
    }
    
    return interactionData.m_fScalingFactor;
}

void ResourceAllocation::notifyResourceConsumed(const std::string& resourceName, double amount)
{
    // Don't record zero or negative consumption
    if (amount <= 0.0)
    {
        return;
    }
    
    // Create or get the resource record
    auto& resource = m_resources[resourceName];
    
    // If this is a dry run
    if (m_curRealRunId == 0)
    {
        // Track consumption during dry run
        resource.m_fConsumed += amount;
        
        // Add this resource to the current interaction's list of consumed resources
        auto currentInteractionIt = m_interactions.begin();
        while (currentInteractionIt != m_interactions.end() && currentInteractionIt->second.m_dryRunId != m_curDryRunId)
        {
            ++currentInteractionIt;
        }
        
        if (currentInteractionIt != m_interactions.end())
        {
            // Add to consumed resources if not already there
            auto& consumedResources = currentInteractionIt->second.m_consumedResourceNames;
            if (std::find(consumedResources.begin(), consumedResources.end(), resourceName) == consumedResources.end())
            {
                consumedResources.push_back(resourceName);
            }
        }
    }
}

void ResourceAllocation::notifyNewRealRun()
{
    assert(m_curRealRunId < m_curDryRunId);
    m_curRealRunId = m_curDryRunId;
}

void ResourceAllocation::updateAvailableResources(const GridCell &cell)
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
