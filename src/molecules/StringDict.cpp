#include "StringDict.h"

// Define static members
std::vector<std::string> StringDict::m_idToString;
std::unordered_map<std::string, StringDict::ID> StringDict::m_stringToId;

void StringDict::initialize()
{
    if (m_idToString.size() > 0)
        return; // already initialized
    
    // Resize vector to accommodate all enum values
    m_idToString.resize(static_cast<size_t>(ID::ORGANELLE_END) + 1);
    
    // Binding surface types
    m_idToString[static_cast<size_t>(ID::eUNKNOWN)] = "UNKNOWN";
    m_idToString[static_cast<size_t>(ID::BS_MEMBRANE)] = "MEMBRANE";
    m_idToString[static_cast<size_t>(ID::BS_CORTEX)] = "CORTEX";
    m_idToString[static_cast<size_t>(ID::BS_CENTROSOME)] = "CENTROSOME";
    
    // PAR proteins
    m_idToString[static_cast<size_t>(ID::PAR_1)] = "PAR-1";
    m_idToString[static_cast<size_t>(ID::PAR_2)] = "PAR-2";
    m_idToString[static_cast<size_t>(ID::PAR_3)] = "PAR-3";
    m_idToString[static_cast<size_t>(ID::PAR_6)] = "PAR-6";
    m_idToString[static_cast<size_t>(ID::PKC_3)] = "PKC-3";
    
    // Cell cycle proteins
    m_idToString[static_cast<size_t>(ID::CDK_1)] = "CDK-1";
    m_idToString[static_cast<size_t>(ID::CDK_2)] = "CDK-2";
    m_idToString[static_cast<size_t>(ID::CYB_1)] = "CYB-1";
    m_idToString[static_cast<size_t>(ID::CCE_1)] = "CCE-1";
    m_idToString[static_cast<size_t>(ID::PLK_1)] = "PLK-1";
    m_idToString[static_cast<size_t>(ID::PLK_4)] = "PLK-4";
    
    // Centrosome proteins
    m_idToString[static_cast<size_t>(ID::GAMMA_TUBULIN)] = "γ-TUBULIN";
    m_idToString[static_cast<size_t>(ID::PERICENTRIN)] = "PERICENTRIN";
    m_idToString[static_cast<size_t>(ID::NINEIN)] = "NINEIN";
    
    // Cell fate specification genes
    m_idToString[static_cast<size_t>(ID::MEX_3)] = "mex-3";
    m_idToString[static_cast<size_t>(ID::SKN_1)] = "skn-1";
    m_idToString[static_cast<size_t>(ID::PAL_1)] = "pal-1";
    m_idToString[static_cast<size_t>(ID::PIE_1)] = "pie-1";
    
    // Organelle types
    m_idToString[static_cast<size_t>(ID::ORGANELLE_NUCLEUS)] = "NUCLEUS";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_MITOCHONDRION)] = "MITOCHONDRION";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_ENDOPLASMIC_RETICULUM)] = "ENDOPLASMIC_RETICULUM";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_SPINDLE)] = "SPINDLE";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_CENTROSOME)] = "CENTROSOME";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_CORTEX)] = "CORTEX_ORGANELLE";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_END)] = "ORGANELLE_END";
    
    // Build reverse mapping
    for (size_t i = 0; i < m_idToString.size(); ++i) {
        if (!m_idToString[i].empty()) {  // Only add non-empty strings
            m_stringToId[m_idToString[i]] = static_cast<ID>(i);
        }
    }
}
