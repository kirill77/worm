#include "pch.h"
#include "GeneWiki.h"
#include "StringDict.h"
#include <stdexcept>

GeneWiki::GeneWiki()
{
    initializeDefaultSequences();
}

GeneWiki& GeneWiki::getInstance()
{
    static GeneWiki instance;
    return instance;
}

const std::string& GeneWiki::getSequence(const std::string& geneName) const
{
    auto it = m_sequences.find(geneName);
    if (it == m_sequences.end())
    {
        throw std::runtime_error("Gene sequence not found: " + geneName);
    }
    return it->second;
}

void GeneWiki::initializeDefaultSequences()
{
    // Cell fate specification genes
    m_sequences[StringDict::idToString(StringDict::ID::PIE_1)] = "ATGCCGAATTCGTCGAATCCG";  // Germline specification
    m_sequences[StringDict::idToString(StringDict::ID::PAL_1)] = "ATGAATTCGCCGAATCCGTCG";  // Posterior fate
    m_sequences[StringDict::idToString(StringDict::ID::SKN_1)] = "ATGCCGTCGAATTCGAATCCG";  // Endoderm specification
    m_sequences[StringDict::idToString(StringDict::ID::MEX_3)] = "ATGTCGCCGAATTCGAATCCG";  // Anterior fate
    
    // Cell division and timing genes
    m_sequences[StringDict::idToString(StringDict::ID::CDK_1)] = "ATGCCGAATTCGTCGAATCCG";  // Cell cycle control
    m_sequences[StringDict::idToString(StringDict::ID::CDK_2)] = "ATGCCGAAGTCGTCGAATCCG";  // CDK-2 transcriptional regulator
    m_sequences[StringDict::idToString(StringDict::ID::CYB_1)] = "ATGAATTCGCCGTCGAATCCG";  // Cyclin B
    m_sequences[StringDict::idToString(StringDict::ID::CCE_1)] = "ATGAAGTTCGCCGAATCCGTC";  // Cyclin E transcriptional regulator
    m_sequences[StringDict::idToString(StringDict::ID::PLK_1)] = "ATGCCGTCGAATTCGAATCCG";  // Polo-like kinase
    
    // Centrosome proteins
    m_sequences[StringDict::idToString(StringDict::ID::GAMMA_TUBULIN)] = "ATGGCCGTCGAATTCCTGACC";  // γ-tubulin
    
    // tRNA genes - these represent the tRNA molecules themselves
    // Start codon tRNA (essential for translation initiation)
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_MET_ATG)] = "ATGGCCAAGCTGAAGTAG";  // Met-ATG initiator tRNA
    
    // Common amino acid tRNAs (high abundance)
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_GLY_GGA)] = "GGATCCAAGCTGGAGTAG";  // Gly-GGA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_GLY_GGT)] = "GGTACCAAGCTGGAGTAG";  // Gly-GGT
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ALA_GCA)] = "GCAAAGCTGAAGTAG";     // Ala-GCA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ALA_GCC)] = "GCCAAGCTGAAGTAG";     // Ala-GCC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_LEU_CTG)] = "CTGGCCAAGCTGAAGTAG";  // Leu-CTG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_LEU_CTC)] = "CTCGCCAAGCTGAAGTAG";  // Leu-CTC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_SER_TCA)] = "TCAAAGCTGAAGTAG";     // Ser-TCA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_SER_TCG)] = "TCGAAGCTGAAGTAG";     // Ser-TCG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_VAL_GTG)] = "GTGGCCAAGCTGAAGTAG";  // Val-GTG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_VAL_GTC)] = "GTCGCCAAGCTGAAGTAG";  // Val-GTC
    
    // Additional essential amino acid tRNAs
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_PRO_CCA)] = "CCAAAGCTGAAGTAG";     // Pro-CCA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_THR_ACA)] = "ACAAAGCTGAAGTAG";     // Thr-ACA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ASP_GAC)] = "GACAAGCTGAAGTAG";     // Asp-GAC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_GLU_GAG)] = "GAGGCCAAGCTGAAGTAG";  // Glu-GAG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_LYS_AAG)] = "AAGGCCAAGCTGAAGTAG";  // Lys-AAG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ARG_CGA)] = "CGAAAGCTGAAGTAG";     // Arg-CGA
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_HIS_CAC)] = "CACAAGCTGAAGTAG";     // His-CAC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_PHE_TTC)] = "TTCAAGCTGAAGTAG";     // Phe-TTC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_TYR_TAC)] = "TACAAGCTGAAGTAG";     // Tyr-TAC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_CYS_TGC)] = "TGCAAGCTGAAGTAG";     // Cys-TGC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_TRP_TGG)] = "TGGAAGCTGAAGTAG";     // Trp-TGG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ASN_AAC)] = "AACAAGCTGAAGTAG";     // Asn-AAC
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_GLN_CAG)] = "CAGAAGCTGAAGTAG";     // Gln-CAG
    m_sequences[StringDict::idToString(StringDict::ID::TRNA_ILE_ATC)] = "ATCAAGCTGAAGTAG";     // Ile-ATC
    
    // Note: These are simplified representative sequences for tRNA genes
    // In reality, tRNA genes have complex secondary structures and processing
    // The sequences here allow for codon matching during translation simulation
} 