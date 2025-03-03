#pragma once

#include "Organelle.h"
#include "DNA.h"
#include "TRNA.h"

class Nucleus : public Organelle
{
private:
    DNA m_genome; // some representation of DNA/genes
    std::vector<TRNA> m_tRNA; // Stores available tRNAs

public:
    void update(double dt) override
    {
        // Simulate transcription process
        // e.g., produce mRNA based on transcriptionRate
    }
};

