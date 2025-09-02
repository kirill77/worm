#pragma once

#include <unordered_map>
#include <string>
#include <assert.h>
#include <vector>

struct StringDict
{
    enum class ID {
        eUNKNOWN,

        // PAR proteins (polarity establishment)
        PAR_1,
        PAR_2,
        PAR_3,
        PAR_6,
        PKC_3,
        
        // Cell cycle proteins
        CDK_1,
        CDK_2,
        CYB_1,
        CCE_1,
        PLK_1,
        PLK_4,
        
        // Centrosome proteins
        GAMMA_TUBULIN,
        PERICENTRIN,
        NINEIN,
        
        // Nucleotides
        ATP,
        
        // tRNA genes (essential set for translation)
        // Start codon
        TRNA_MET_ATG,        // Methionine initiator tRNA
        
        // Common amino acids (high abundance needed)
        TRNA_GLY_GGA,        // Glycine (preferred codon in C. elegans)
        TRNA_GLY_GGT,        // Glycine (second choice)
        TRNA_ALA_GCA,        // Alanine (preferred codon)
        TRNA_ALA_GCC,        // Alanine (second choice)
        TRNA_LEU_CTG,        // Leucine (highly preferred in C. elegans)
        TRNA_LEU_CTC,        // Leucine (second choice)
        TRNA_SER_TCA,        // Serine (common codon)
        TRNA_SER_TCG,        // Serine (alternative)
        TRNA_VAL_GTG,        // Valine (preferred)
        TRNA_VAL_GTC,        // Valine (alternative)
        
        // Less common but essential amino acids
        TRNA_PRO_CCA,        // Proline
        TRNA_THR_ACA,        // Threonine  
        TRNA_ASP_GAC,        // Aspartic acid
        TRNA_GLU_GAG,        // Glutamic acid
        TRNA_LYS_AAG,        // Lysine
        TRNA_ARG_CGA,        // Arginine
        TRNA_HIS_CAC,        // Histidine
        TRNA_PHE_TTC,        // Phenylalanine
        TRNA_TYR_TAC,        // Tyrosine
        TRNA_CYS_TGC,        // Cysteine
        TRNA_TRP_TGG,        // Tryptophan (single codon)
        TRNA_ASN_AAC,        // Asparagine
        TRNA_GLN_CAG,        // Glutamine  
        TRNA_ILE_ATC,        // Isoleucine
        
        // Cell fate specification genes
        MEX_3,
        SKN_1,
        PAL_1,
        PIE_1,
        
        // Organelle types (must be contiguous for vector indexing)
        ORGANELLE_START,
        ORGANELLE_NUCLEUS = ORGANELLE_START,
        ORGANELLE_MITOCHONDRION,
        ORGANELLE_ENDOPLASMIC_RETICULUM,
        ORGANELLE_SPINDLE,
        ORGANELLE_CENTROSOME,
        ORGANELLE_CORTEX,
        ORGANELLE_END,
        
        // TODO: add all IDs here
    };

    static void initialize();
    static const std::string& idToString(ID id)
    {
        assert(m_idToString.size() > 0);
        return m_idToString[static_cast<size_t>(id)];
    }
    static ID stringToId(const std::string& s);

private:
    static std::vector<std::string> m_idToString;
    static std::unordered_map<std::string, ID> m_stringToId;
};

