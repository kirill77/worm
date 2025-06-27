#include "StringDict.h"

// Define static members
std::vector<std::string> StringDict::m_idToString;
std::unordered_map<std::string, StringDict::ID> StringDict::m_stringToId;

void StringDict::initialize()
{
    if (m_idToString.size() > 0)
        return; // already initialized
    
    // Add binding surface strings
    m_idToString.push_back("UNKNOWN");  // eUNKNOWN
    m_idToString.push_back("MEMBRANE"); // BS_MEMBRANE
    m_idToString.push_back("CORTEX");   // BS_CORTEX
    m_idToString.push_back("CENTROSOME"); // BS_CENTROSOME
    
    // Add PAR proteins
    m_idToString.push_back("PAR-1");    // PAR_1
    m_idToString.push_back("PAR-2");    // PAR_2
    m_idToString.push_back("PAR-3");    // PAR_3
    m_idToString.push_back("PAR-6");    // PAR_6
    m_idToString.push_back("PKC-3");    // PKC_3
    
    // Add cell cycle proteins
    m_idToString.push_back("CDK-1");    // CDK_1
    m_idToString.push_back("CDK-2");    // CDK_2
    m_idToString.push_back("CYB-1");    // CYB_1
    m_idToString.push_back("CCE-1");    // CCE_1
    m_idToString.push_back("PLK-1");    // PLK_1
    m_idToString.push_back("PLK-4");    // PLK_4
    
    // Add centrosome proteins
    m_idToString.push_back("γ-TUBULIN"); // GAMMA_TUBULIN
    m_idToString.push_back("PERICENTRIN"); // PERICENTRIN
    m_idToString.push_back("NINEIN");   // NINEIN
    
    // Add cell fate specification genes
    m_idToString.push_back("mex-3");    // MEX_3
    m_idToString.push_back("skn-1");    // SKN_1
    m_idToString.push_back("pal-1");    // PAL_1
    m_idToString.push_back("pie-1");    // PIE_1
    
    // Add organelle types (ORGANELLE_START is same as ORGANELLE_NUCLEUS)
    m_idToString.push_back("NUCLEUS");               // ORGANELLE_NUCLEUS/ORGANELLE_START
    m_idToString.push_back("MITOCHONDRION");         // ORGANELLE_MITOCHONDRION
    m_idToString.push_back("ENDOPLASMIC_RETICULUM"); // ORGANELLE_ENDOPLASMIC_RETICULUM
    m_idToString.push_back("SPINDLE");               // ORGANELLE_SPINDLE
    m_idToString.push_back("CENTROSOME");            // ORGANELLE_CENTROSOME
    m_idToString.push_back("ORGANELLE_END");         // ORGANELLE_END
    
    // Build reverse mapping
    for (size_t i = 0; i < m_idToString.size(); ++i) {
        m_stringToId[m_idToString[i]] = static_cast<ID>(i);
    }
}
