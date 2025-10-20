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

    // Microtubule dynamic instability (catastrophe-related parameters)
    // Note: These affect catastrophe probability directly or via cap/tip modulators
    constexpr double MT_CATASTROPHE_BASE_FREE            = 0.12;  // s^-1 baseline (free tip)
    constexpr double MT_CATASTROPHE_BASE_CONTACT         = 0.60;  // s^-1 when contacting cortex
    constexpr double MT_K_TUB                            = 50.0;  // molecules/µm^3 (tubulin saturation scale)
    constexpr double MT_K_AIR                            = 10.0;  // molecules/µm^3 (AIR-1 modulation scale)
    constexpr double MT_CAP_HYDROLYSIS_RATE_S            = 0.03;  // s^-1 (used to drain cap length proxy)
    constexpr double MT_CAP_DEPLETED_THRESHOLD_MICROM    = 0.02;  // µm (cap length threshold for high catastrophe)
    constexpr double MT_CAP_DEPLETION_CATASTROPHE_MULT   = 3.0;   // catastrophe multiplier when cap depleted

    // Microtubule growth/shrink speeds
    constexpr double MT_VGROW_MAX_UM_PER_S = 0.45;  // µm/s at high tubulin
    constexpr double MT_VSHRINK_UM_PER_S   = 0.09;  // µm/s
    
    // Microtubule segmentation (for bendable microtubules)
    constexpr double MT_SEGMENT_LENGTH_MICROM = 0.25;  // µm per segment (250 nm)
    
    // Cortical binding (dynein-mediated MT-cortex anchoring)
    constexpr double MT_CORTEX_K_BIND      = 100.0;  // molecules/µm^3 (dynein concentration for half-max binding)
    constexpr double MT_CORTEX_BIND_RATE   = 2.0;    // s^-1 (attempt rate when in contact)
    constexpr double MT_CORTEX_UNBIND_BASE = 0.1;    // s^-1 (base unbinding rate)
    constexpr double MT_CORTEX_UNBIND_WEAK = 1.0;    // s^-1 (unbinding rate for weak binding)
    
    // Dynein pulling force (physics)
    constexpr double DYNEIN_PULLING_FORCE_PICONEWTONS = 5.0;  // pN (picoNewtons per bound microtubule)
}


