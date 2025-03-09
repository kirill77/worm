#pragma once

#include "Organelle.h"
#include "math/vector.h"
#include "CellTypes.h"

class Spindle : public Organelle
{
private:
    float3 m_position;        // Center position of the spindle
    float3 m_orientation;     // Normalized direction vector (from - to + end)
    float m_length;          // Current length of the spindle
    bool m_isAssembled;       // Whether the spindle is fully assembled
    CellType m_cellType;      // Type of cell this spindle is in
    
    // Microtubule parameters
    static constexpr float INITIAL_LENGTH = 0.2f;       // Initial spindle length
    static constexpr float MAX_LENGTH = 0.8f;           // Maximum spindle length
    static constexpr float GROWTH_RATE = 0.1f;          // Growth rate during assembly
    static constexpr float CORTICAL_FORCE = 0.5f;       // Force from cortical pulling
    static constexpr float ROTATION_RATE = 0.2f;        // Rate of orientation change
    static constexpr float ATP_PER_GROWTH = 5.0f;       // ATP cost for length increase
    static constexpr float ATP_PER_ROTATION = 2.0f;     // ATP cost for orientation change

public:
    Spindle(CellType type = CellType::Zygote)
        : m_position(0.0f, 0.0f, 0.0f)
        , m_orientation(1.0f, 0.0f, 0.0f)  // Start aligned with x-axis
        , m_length(INITIAL_LENGTH)
        , m_isAssembled(false)
        , m_cellType(type)
    {}

    void update(double dt, Cell& cell, Medium& medium) override;

    // Getters for spindle state
    bool isAssembled() const { return m_isAssembled; }
    float3 getPosition() const { return m_position; }
    float3 getOrientation() const { return m_orientation; }
    float getLength() const { return m_length; }
    CellType getCellType() const { return m_cellType; }
    
    // Get positions of spindle poles
    float3 getMinusPole() const { return m_position - m_orientation * (m_length * 0.5f); }
    float3 getPlusPole() const { return m_position + m_orientation * (m_length * 0.5f); }

private:
    // Helper functions for spindle dynamics
    void updateAssembly(double dt, Cell& cell, Medium& medium);
    void updatePosition(double dt, Cell& cell, Medium& medium);
    void updateOrientation(double dt, Cell& cell, Medium& medium);
    
    // Calculate forces based on cell type
    float3 calculateCorticalForces(const Medium& medium) const;
    float3 calculatePARBasedForces(const Medium& medium) const;  // For P lineage
    float3 calculateDefaultForces(const Medium& medium) const;   // For somatic cells
    
    // Calculate preferred orientation based on cell type
    float3 calculatePreferredOrientation(const Medium& medium) const;
}; 