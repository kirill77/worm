#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <assert.h>
#include "Molecule.h"

// There is a list of resources and a list of interactions that use those resources. If
// we simply apply the interactions - the interactions that are being applied first will have
// an unfair advantage because they can consume all resources they want. This class is used
// to fix that. It assumes the interactions will be applied in two passes. First pass (dry run)
// gathers the information about used resources, and second pass distributes the resources
// fairly between all interactions
class ResourceDistributor
{
public:
    ResourceDistributor();
    ~ResourceDistributor();

    void notifyNewDryRun(const class GridCell& cell);

    // if this returns false - you can skip this interaction
    bool notifyNewInteractionStarting(const class MoleculeInteraction &interaction);

    double getAvailableResource(const Molecule& molecule);

    void notifyResourceWanted(const Molecule& molecule, double m_fNumber);

    void notifyNewRealRun();

    bool isDryRun() const { return m_curDryRunId > m_curRealRunId; }

private:
    void updateAvailableResources(const GridCell& cell);

    uint64_t m_curDryRunId = 0, m_curRealRunId = 0;

    struct ResourceData
    {
        uint64_t m_dryRunId = 0;
        double m_fRequested = 0, m_fAvailable = 0;
        double computeScalingFactor()
        {
            assert(m_fRequested >= 0 && m_fAvailable >= 0);
            return m_fAvailable >= m_fRequested ? 1 : m_fAvailable / m_fRequested;
        }
    };
    std::unordered_map<Molecule, ResourceData> m_resources;

    struct InteractionData
    {
        uint64_t m_lastValidDryRunId = 0;
        double m_fScalingFactor = 1.0;
        std::vector<Molecule> m_requestedMolecules;
    };
    std::unordered_map<const MoleculeInteraction*, InteractionData> m_interactions;

    InteractionData* m_pCurInteraction = nullptr;
}; 