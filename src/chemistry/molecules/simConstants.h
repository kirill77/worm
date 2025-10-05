#pragma once

namespace MoleculeConstants
{
    // Pol II transcription regulation (DNA::updateTranscriptionalRegulation)
    constexpr double TRANSCRIPTION_BASAL_RATE = 0.05;         // Basal transcription rate
    constexpr double TRANSCRIPTION_MAX_ACTIVATED_RATE = 0.8;   // Max activated transcription rate
    constexpr double TF_ACTIVITY_K = 250000.0;                 // Denominator in TF activity term

    // Pol III tRNA production (DNA::transcribeAll)
    constexpr double TRNA_POLIII_PRODUCTION_MULTIPLIER = 1000000.0; // Current diagnostic multiplier

    // Nuclear transport and envelope thresholds (Nucleus)
    constexpr double NUCLEAR_TF_IMPORT_RATE = 0.1;             // Fraction per step imported to nucleus
    constexpr double ENVELOPE_TRANSCRIBE_THRESHOLD = 0.8;      // Transcription allowed if envelope >= this
    constexpr double ENVELOPE_EXPORT_THRESHOLD = 0.5;          // Export allowed if envelope >= this

    // Tubulin gene expression defaults (used in Worm gene initialization)
    constexpr double ALPHA_TUBULIN_EXPRESSION_RATE = 50000.0;
    constexpr double BETA_TUBULIN_EXPRESSION_RATE  = 50000.0;
    constexpr double ALPHA_TUBULIN_BASAL_LEVEL     = 0.2;
    constexpr double BETA_TUBULIN_BASAL_LEVEL      = 0.2;
}


