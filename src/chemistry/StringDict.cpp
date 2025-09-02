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
    
    // PAR proteins (polarity establishment)
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
    
    // Nucleotides
    m_idToString[static_cast<size_t>(ID::ATP)] = "ATP";
    
    // tRNA genes (essential set for translation)
    // Start codon
    m_idToString[static_cast<size_t>(ID::TRNA_MET_ATG)] = "tRNA-Met-ATG";
    
    // Common amino acids (high abundance)
    m_idToString[static_cast<size_t>(ID::TRNA_GLY_GGA)] = "tRNA-Gly-GGA";
    m_idToString[static_cast<size_t>(ID::TRNA_GLY_GGT)] = "tRNA-Gly-GGT";
    m_idToString[static_cast<size_t>(ID::TRNA_ALA_GCA)] = "tRNA-Ala-GCA";
    m_idToString[static_cast<size_t>(ID::TRNA_ALA_GCC)] = "tRNA-Ala-GCC";
    m_idToString[static_cast<size_t>(ID::TRNA_LEU_CTG)] = "tRNA-Leu-CTG";
    m_idToString[static_cast<size_t>(ID::TRNA_LEU_CTC)] = "tRNA-Leu-CTC";
    m_idToString[static_cast<size_t>(ID::TRNA_SER_TCA)] = "tRNA-Ser-TCA";
    m_idToString[static_cast<size_t>(ID::TRNA_SER_TCG)] = "tRNA-Ser-TCG";
    m_idToString[static_cast<size_t>(ID::TRNA_VAL_GTG)] = "tRNA-Val-GTG";
    m_idToString[static_cast<size_t>(ID::TRNA_VAL_GTC)] = "tRNA-Val-GTC";
    
    // Less common but essential amino acids  
    m_idToString[static_cast<size_t>(ID::TRNA_PRO_CCA)] = "tRNA-Pro-CCA";
    m_idToString[static_cast<size_t>(ID::TRNA_THR_ACA)] = "tRNA-Thr-ACA";
    m_idToString[static_cast<size_t>(ID::TRNA_ASP_GAC)] = "tRNA-Asp-GAC";
    m_idToString[static_cast<size_t>(ID::TRNA_GLU_GAG)] = "tRNA-Glu-GAG";
    m_idToString[static_cast<size_t>(ID::TRNA_LYS_AAG)] = "tRNA-Lys-AAG";
    m_idToString[static_cast<size_t>(ID::TRNA_ARG_CGA)] = "tRNA-Arg-CGA";
    m_idToString[static_cast<size_t>(ID::TRNA_HIS_CAC)] = "tRNA-His-CAC";
    m_idToString[static_cast<size_t>(ID::TRNA_PHE_TTC)] = "tRNA-Phe-TTC";
    m_idToString[static_cast<size_t>(ID::TRNA_TYR_TAC)] = "tRNA-Tyr-TAC";
    m_idToString[static_cast<size_t>(ID::TRNA_CYS_TGC)] = "tRNA-Cys-TGC";
    m_idToString[static_cast<size_t>(ID::TRNA_TRP_TGG)] = "tRNA-Trp-TGG";
    m_idToString[static_cast<size_t>(ID::TRNA_ASN_AAC)] = "tRNA-Asn-AAC";
    m_idToString[static_cast<size_t>(ID::TRNA_GLN_CAG)] = "tRNA-Gln-CAG";
    m_idToString[static_cast<size_t>(ID::TRNA_ILE_ATC)] = "tRNA-Ile-ATC";
    
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
    m_idToString[static_cast<size_t>(ID::ORGANELLE_CORTEX)] = "CORTEX";
    m_idToString[static_cast<size_t>(ID::ORGANELLE_END)] = "ORGANELLE_END";
    
    // Build reverse mapping
    for (size_t i = 0; i < m_idToString.size(); ++i) {
        if (!m_idToString[i].empty()) {  // Only add non-empty strings
            m_stringToId[m_idToString[i]] = static_cast<ID>(i);
        }
    }
}

StringDict::ID StringDict::stringToId(const std::string& s)
{
    assert(m_stringToId.size() > 0);
    auto it = m_stringToId.find(s);
    assert(it != m_stringToId.end() && "String not found in StringDict - make sure to call initialize() and add the string");
    return (it == m_stringToId.end()) ? ID::eUNKNOWN : it->second;
}