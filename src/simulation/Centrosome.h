#pragma once

#include "Organelle.h"
#include "geometry/vectors/vector.h"

class Centrosome : public Organelle
{
private:
    float3 m_position;  // Position of the centrosome in the cell
    bool m_isDuplicated;  // Whether the centrosome has duplicated
    double m_duplicationTime;  // Time when duplication occurred

public:
    /**
     * Constructor for Centrosome
     * 
     * @param pCell Weak pointer to the cell containing this centrosome
     * @param position Initial position of the centrosome
     */
    Centrosome(std::weak_ptr<Cell> pCell, const float3& position = float3(0, 0, 0));
    
    /**
     * Update function for centrosome behavior
     * 
     * @param dt Time step in seconds
     * @param cell Reference to the cell containing this centrosome
     */
    void update(double dt, Cell& cell) override;
    
    /**
     * Get the position of the centrosome
     * 
     * @return Position vector
     */
    const float3& getPosition() const { return m_position; }
    
    /**
     * Set the position of the centrosome
     * 
     * @param position New position
     */
    void setPosition(const float3& position) { m_position = position; }
    
    /**
     * Check if the centrosome has duplicated
     * 
     * @return True if duplicated, false otherwise
     */
    bool isDuplicated() const { return m_isDuplicated; }
    
    /**
     * Trigger centrosome duplication
     */
    void duplicate();
};

