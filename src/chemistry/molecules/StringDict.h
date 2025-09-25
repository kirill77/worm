#pragma once

#include <unordered_map>
#include <string>
#include <assert.h>
#include <vector>

struct StringDict
{
    enum class ID {
        eUNKNOWN,

        // Define contiguous gene range for scanning genes/mRNAs
        GENES_START,

        // PAR proteins (polarity establishment)
        PAR_1 = GENES_START,
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

        // Centrosome structural proteins (genes)
        GAMMA_TUBULIN,
        PERICENTRIN,
        NINEIN,

        GENES_END,

        // Non-gene molecules
        // Nucleotides
        ATP,
        GTP,
        GDP,

        // Endoplasmic reticulum molecules
        ER_PROTEIN,        // ER-synthesized proteins
        ER_LIPID,          // ER-synthesized lipids

        // Rho/CDC42 pathway and cortex mechanics (proteins and proxies)
        RHO_1,             // Total RHO-1
        RHO_1_GTP,         // Active RHO-1
        RHO_1_GDP,         // Inactive RHO-1
        ECT_2,             // RhoGEF ECT-2
        CHIN_1,            // RhoGAP (CHIN-1 proxy)
        CDC_42,            // Total CDC-42
        CDC_42_GTP,        // Active CDC-42
        CDC_42_GDP,        // Inactive CDC-42
        NMY_2,             // Non-muscle myosin II (contractility proxy)
        AIR_1,             // Aurora A kinase at centrosome/MTs
        MT_PLUS_END_CONTACT, // Marker for MT plus-end cortex contact (proxy species)
        SPD_2,             // PCM scaffold CEP192 homolog
        SPD_5,             // PCM scaffold
        PCMD_1,            // Centriolar–PCM bridge protein
        ZYG_1,             // Plk4-like duplication kinase
        SAS_6,             // Cartwheel protein
        SAS_5,             // STIL/Ana2-like
        SAS_4,             // Centriole outer layer assembly

        // Phosphorylated PAR proteins
        PAR_1_P,          // Phosphorylated PAR-1
        PAR_2_P,          // Phosphorylated PAR-2
        PAR_3_P,          // Phosphorylated PAR-3

        // Protein complexes
        PAR_3_PAR_6,      // PAR-3:PAR-6 complex
        PAR_6_PKC_3,      // PAR-6:PKC-3 complex
        PAR_1_CORTEX,     // PAR-1:CORTEX complex
        PAR_2_CORTEX,     // PAR-2:CORTEX complex
        PAR_3_CORTEX,     // PAR-3:CORTEX complex

        // Charged tRNA variants (for distinguishing charged vs uncharged)
        // Start codon - charged
        TRNA_MET_ATG_CHARGED,
        
        // Common amino acids - charged
        TRNA_GLY_GGA_CHARGED,
        TRNA_GLY_GGT_CHARGED,
        TRNA_ALA_GCA_CHARGED,
        TRNA_ALA_GCC_CHARGED,
        TRNA_LEU_CTG_CHARGED,
        TRNA_LEU_CTC_CHARGED,
        TRNA_SER_TCA_CHARGED,
        TRNA_SER_TCG_CHARGED,
        TRNA_VAL_GTG_CHARGED,
        TRNA_VAL_GTC_CHARGED,
        
        // Less common but essential amino acids - charged
        TRNA_PRO_CCA_CHARGED,
        TRNA_THR_ACA_CHARGED,
        TRNA_ASP_GAC_CHARGED,
        TRNA_GLU_GAG_CHARGED,
        TRNA_LYS_AAG_CHARGED,
        TRNA_ARG_CGA_CHARGED,
        TRNA_HIS_CAC_CHARGED,
        TRNA_PHE_TTC_CHARGED,
        TRNA_TYR_TAC_CHARGED,
        TRNA_CYS_TGC_CHARGED,
        TRNA_TRP_TGG_CHARGED,
        TRNA_ASN_AAC_CHARGED,
        TRNA_GLN_CAG_CHARGED,
        TRNA_ILE_ATC_CHARGED,

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

