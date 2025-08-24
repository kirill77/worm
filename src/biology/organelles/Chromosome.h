#pragma once

#include "geometry/vectors/vector.h"
#include "chemistry/DNA.h"
#include <memory>

class Chromosome
{
private:
    std::shared_ptr<DNA> m_pDNA;         // Genetic material in this chromosome
    float m_fCondensation;               // 0.0 = relaxed, 1.0 = fully condensed
    float3 m_position;                   // Position relative to cell center
    bool m_bIsAttached;                  // Whether attached to spindle microtubules
    bool m_bIsSeparated;                 // Whether sister chromatids have separated
    float3 m_attachmentPoint;            // Point where spindle microtubules attach (kinetochore)

    static constexpr float fCONDENSATION_RATE = 0.2f;     // Rate of chromosome condensation
    static constexpr float fDECONDENSATION_RATE = 0.3f;   // Rate of chromosome decondensation
    static constexpr float fSEPARATION_DISTANCE = 0.1f;   // Distance between separated chromatids

public:
    Chromosome(std::shared_ptr<DNA> pDNA = nullptr)
        : m_pDNA(pDNA)
        , m_fCondensation(0.0f)
        , m_position(0.0f, 0.0f, 0.0f)
        , m_bIsAttached(false)
        , m_bIsSeparated(false)
        , m_attachmentPoint(0.0f, 0.0f, 0.0f)
    {}

    // Main update function
    void update(double fDt, class Cell& cell, class Medium& medium);

    // Chromosome state changes
    void condense(float fDt);      // Called during prophase
    void decondense(float fDt);    // Called during telophase
    void separate();               // Called during anaphase
    
    // Spindle interaction
    bool tryAttachToSpindle(const class Spindle& spindle);
    void moveAlongSpindle(const class Spindle& spindle, float fDt);
    
    // DNA-related functions
    void setDNA(std::shared_ptr<DNA> pDNA) { m_pDNA = pDNA; }
    std::shared_ptr<DNA> getDNA() const { return m_pDNA; }
    std::vector<std::shared_ptr<MRNA>> transcribe(double fDt, const class GridCell& nuclearCompartment) const;
    
    // Getters
    float getCondensation() const { return m_fCondensation; }
    float3 getPosition() const { return m_position; }
    bool isAttached() const { return m_bIsAttached; }
    bool isSeparated() const { return m_bIsSeparated; }
    bool isFullyCondensed() const { return m_fCondensation > 0.95f; }
    bool isFullyDecondensed() const { return m_fCondensation < 0.05f; }
}; 