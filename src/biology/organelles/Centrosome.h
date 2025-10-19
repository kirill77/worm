#pragma once

#include <memory>
#include <vector>
#include "Organelle.h"
#include "physics/PhysCentrosome.h"

struct Y_TuRC;

struct Centrosome : public Organelle
{
private:
    std::shared_ptr<PhysCentrosome> m_pPhysCentrosome;
    bool m_isDuplicated;  // Whether the centrosome has duplicated
    double m_duplicationTime;  // Time when duplication occurred
    float m_fPCMRadiusMicroM;  // PCM (Pericentriolar Material) radius in micrometers
    // Simple PCM maturation proxy in [0,1] to drive γ-tubulin recruitment capacity
    double m_pcmMaturation = 0.1;
    // Bound γ-tubulin concentration proxy at centrosome (molecules/µm^3), grid-agnostic
    double m_gammaBoundConc = 0.0;

public:
    /**
     * Constructor for Centrosome
     * 
     * @param pCell Weak pointer to the cell containing this centrosome
     * @param vNormalizedPos Initial position of the centrosome in normalized coordinates (-1, 1)
     */
    Centrosome(std::weak_ptr<Cell> pCell, const float3& vNormalizedPos = float3(0, 0, 0));
    
    /**
     * Update function for centrosome behavior
     * 
     * @param dt Time step in seconds
     * @param cell Reference to the cell containing this centrosome
     */
    void update(double dt, Cell& cell) override;
    
    /**
     * Get the position of the centrosome in parent (cell) space
     * 
     * @return Position vector in normalized coordinates (-1, 1)
     */
    float3 getNormalizedPosition() const { return m_pPhysCentrosome->getToNormalizedCell().m_translation; }

    /**
     * Check if the centrosome has duplicated
     * 
     * @return True if duplicated, false otherwise
     */
    bool isDuplicated() const { return m_isDuplicated; }

private:
    // Update PCM maturation based on local SPD-2/5 and kinase activity at the centrosome
    void updatePCMMaturation(double dt, const Cell& cell, Medium& internalMedium);
    // Update local gamma-tubulin enrichment and adjust Y_TuRC ring complex count
    void updateGammaAndRingComplexes(double dt, const Cell& cell, Medium& internalMedium);

public:
    /**
     * Get the PCM (Pericentriolar Material) radius
     * 
     * @return PCM radius in micrometers
     */
    float getPCMRadius() const { return m_fPCMRadiusMicroM; }
    
    /**
     * Get the ring complexes associated with this centrosome
     * 
     * @return Vector of Y_TuRC ring complexes
     */
    const std::vector<std::shared_ptr<PhysMicrotubule>>& getRingComplexes() const { return m_pPhysCentrosome->getMicrotubules(); }
};

