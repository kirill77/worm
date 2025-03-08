#pragma once

#include "Organelle.h"
#include "DNA.h"
#include <memory>
#include <vector>
#include <functional>

// Forward declarations
class MRNA;

class Nucleus : public Organelle
{
private:
    std::shared_ptr<DNA> m_pDNA;
    double m_envelopeIntegrity;  // 1.0 = intact, 0.0 = broken down
    static constexpr double ENVELOPE_BREAKDOWN_RATE = 0.2;  // Rate of nuclear envelope breakdown
    static constexpr double ENVELOPE_REFORM_RATE = 0.5;    // Rate of nuclear envelope reformation

public:
    Nucleus(std::shared_ptr<DNA> pDNA)
        : m_pDNA(pDNA), m_envelopeIntegrity(1.0) {}

    // Update function now takes cell cycle state and a callback for adding mRNAs
    void update(double dt, CellCycleState cellState, 
                std::function<void(std::shared_ptr<MRNA>)> addMRNA) override;

    double getEnvelopeIntegrity() const { return m_envelopeIntegrity; }
};

