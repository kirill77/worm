#pragma once

#include "math/vector.h"

class Chromosome
{
private:
    float m_condensation;     // 0.0 = relaxed, 1.0 = fully condensed
    float3 m_position;        // Position relative to cell center
    bool m_isAttached;        // Whether attached to spindle microtubules
    bool m_isSeparated;       // Whether sister chromatids have separated
    float3 m_attachmentPoint; // Point where spindle microtubules attach (kinetochore)

    static constexpr float CONDENSATION_RATE = 0.2f;    // Rate of chromosome condensation
    static constexpr float DECONDENSATION_RATE = 0.3f;  // Rate of chromosome decondensation
    static constexpr float SEPARATION_DISTANCE = 0.1f;  // Distance between separated chromatids

public:
    Chromosome()
        : m_condensation(0.0f)
        , m_position(0.0f, 0.0f, 0.0f)
        , m_isAttached(false)
        , m_isSeparated(false)
        , m_attachmentPoint(0.0f, 0.0f, 0.0f)
    {}

    // Main update function
    void update(double dt, class Cell& cell, class Medium& medium);

    // Chromosome state changes
    void condense(float dt);      // Called during prophase
    void decondense(float dt);    // Called during telophase
    void separate();              // Called during anaphase
    
    // Spindle interaction
    bool tryAttachToSpindle(const class Spindle& spindle);
    void moveAlongSpindle(const class Spindle& spindle, float dt);
    
    // Getters
    float getCondensation() const { return m_condensation; }
    float3 getPosition() const { return m_position; }
    bool isAttached() const { return m_isAttached; }
    bool isSeparated() const { return m_isSeparated; }
    bool isFullyCondensed() const { return m_condensation > 0.95f; }
    bool isFullyDecondensed() const { return m_condensation < 0.05f; }
}; 