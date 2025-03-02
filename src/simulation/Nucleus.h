#pragma once

#include "Organelle.h"
#include "DNA.h"

class Nucleus : public Organelle
{
private:
    DNA m_genome; // some representation of DNA/genes
    double m_fTranscriptionRate;

public:
    void update(double dt) override
    {
        // Simulate transcription process
        // e.g., produce mRNA based on transcriptionRate
    }
};

