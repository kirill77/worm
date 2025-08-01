#pragma once

#include <memory>
#include "Organelle.h"
#include "geometry/vectors/affine.h"

struct Y_TuRC;

struct Centrosome : public Organelle
{
private:
    affine3 m_mToParent;  // Transform from centrosome space to parent (cell) space
    bool m_isDuplicated;  // Whether the centrosome has duplicated
    double m_duplicationTime;  // Time when duplication occurred
    float m_fPCMRadiusMicroM;  // PCM (Pericentriolar Material) radius in micrometers
    std::vector<std::shared_ptr<Y_TuRC>> m_pRingComplexes;

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
    float3 getNormalizedPosition() const { return m_mToParent.m_translation; }

    /**
     * Get the full transform from centrosome space to parent (cell) space
     * 
     * @return Affine transform matrix
     */
    const affine3& getToParentTransform() const { return m_mToParent; }

    /**
     * Set the position of the centrosome in parent (cell) space
     * 
     * @param vNormalizedPos New position in normalized coordinates (-1, 1)
     */
    void setNormalizedPosition(const float3& vNormalizedPos) { m_mToParent.m_translation = vNormalizedPos; }

    /**
     * Set the full transform from centrosome space to parent (cell) space
     * 
     * @param mTransform Affine transform matrix
     */
    void setToParentTransform(const affine3& mTransform) { m_mToParent = mTransform; }

    /**
     * Check if the centrosome has duplicated
     * 
     * @return True if duplicated, false otherwise
     */
    bool isDuplicated() const { return m_isDuplicated; }
    
    /**
     * Get the PCM (Pericentriolar Material) radius
     * 
     * @return PCM radius in micrometers
     */
    float getPCMRadius() const { return m_fPCMRadiusMicroM; }
    
    /**
     * Set the PCM (Pericentriolar Material) radius
     * 
     * @param radius PCM radius in micrometers
     */
    void setPCMRadius(float radius) { m_fPCMRadiusMicroM = radius; }
    
    /**
     * Get the ring complexes associated with this centrosome
     * 
     * @return Vector of Y_TuRC ring complexes
     */
    const std::vector<std::shared_ptr<Y_TuRC>>& getRingComplexes() const { return m_pRingComplexes; }
    
    /**
     * Trigger centrosome duplication
     */
    void duplicate();
};

