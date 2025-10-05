#include "TRNA.h"
#include "utils/log/ILog.h"
#include <cassert>
#include <algorithm>

bool TRNA::isChargedTRNA(StringDict::ID id)
{
    return id >= StringDict::ID::TRNA_MET_ATG_CHARGED && id <= StringDict::ID::TRNA_ILE_ATC_CHARGED;
}

bool TRNA::isTRNAGeneId(StringDict::ID id)
{
    return id >= StringDict::ID::TRNA_GENES_START && id <= StringDict::ID::TRNA_GENES_END;
}

StringDict::ID TRNA::getChargedVariant(StringDict::ID unchargedID)
{
    switch (unchargedID)
    {
        // Start codon
        case StringDict::ID::TRNA_MET_ATG:          return StringDict::ID::TRNA_MET_ATG_CHARGED;
        
        // Common amino acids
        case StringDict::ID::TRNA_GLY_GGA:          return StringDict::ID::TRNA_GLY_GGA_CHARGED;
        case StringDict::ID::TRNA_GLY_GGT:          return StringDict::ID::TRNA_GLY_GGT_CHARGED;
        case StringDict::ID::TRNA_ALA_GCA:          return StringDict::ID::TRNA_ALA_GCA_CHARGED;
        case StringDict::ID::TRNA_ALA_GCC:          return StringDict::ID::TRNA_ALA_GCC_CHARGED;
        case StringDict::ID::TRNA_LEU_CTG:          return StringDict::ID::TRNA_LEU_CTG_CHARGED;
        case StringDict::ID::TRNA_LEU_CTC:          return StringDict::ID::TRNA_LEU_CTC_CHARGED;
        case StringDict::ID::TRNA_SER_TCA:          return StringDict::ID::TRNA_SER_TCA_CHARGED;
        case StringDict::ID::TRNA_SER_TCG:          return StringDict::ID::TRNA_SER_TCG_CHARGED;
        case StringDict::ID::TRNA_VAL_GTG:          return StringDict::ID::TRNA_VAL_GTG_CHARGED;
        case StringDict::ID::TRNA_VAL_GTC:          return StringDict::ID::TRNA_VAL_GTC_CHARGED;
        
        // Less common but essential amino acids
        case StringDict::ID::TRNA_PRO_CCA:          return StringDict::ID::TRNA_PRO_CCA_CHARGED;
        case StringDict::ID::TRNA_THR_ACA:          return StringDict::ID::TRNA_THR_ACA_CHARGED;
        case StringDict::ID::TRNA_ASP_GAC:          return StringDict::ID::TRNA_ASP_GAC_CHARGED;
        case StringDict::ID::TRNA_GLU_GAG:          return StringDict::ID::TRNA_GLU_GAG_CHARGED;
        case StringDict::ID::TRNA_LYS_AAG:          return StringDict::ID::TRNA_LYS_AAG_CHARGED;
        case StringDict::ID::TRNA_ARG_CGA:          return StringDict::ID::TRNA_ARG_CGA_CHARGED;
        case StringDict::ID::TRNA_HIS_CAC:          return StringDict::ID::TRNA_HIS_CAC_CHARGED;
        case StringDict::ID::TRNA_PHE_TTC:          return StringDict::ID::TRNA_PHE_TTC_CHARGED;
        case StringDict::ID::TRNA_TYR_TAC:          return StringDict::ID::TRNA_TYR_TAC_CHARGED;
        case StringDict::ID::TRNA_CYS_TGC:          return StringDict::ID::TRNA_CYS_TGC_CHARGED;
        case StringDict::ID::TRNA_TRP_TGG:          return StringDict::ID::TRNA_TRP_TGG_CHARGED;
        case StringDict::ID::TRNA_ASN_AAC:          return StringDict::ID::TRNA_ASN_AAC_CHARGED;
        case StringDict::ID::TRNA_GLN_CAG:          return StringDict::ID::TRNA_GLN_CAG_CHARGED;
        case StringDict::ID::TRNA_ILE_ATC:          return StringDict::ID::TRNA_ILE_ATC_CHARGED;
        
        default:
            const std::string& name = StringDict::idToString(unchargedID);
            LOG_ERROR("getChargedVariant called with invalid ID: %s (ID %d) is not an uncharged tRNA", 
                     name.c_str(), static_cast<int>(unchargedID));
            assert(false && "getChargedVariant called with non-uncharged-tRNA ID");
            return StringDict::ID::eUNKNOWN;
    }
}

// Returns the anticodon sequence for a given tRNA molecule
// 
// BIOLOGICAL BACKGROUND:
// Each tRNA molecule has a specific 3-nucleotide anticodon sequence in its anticodon loop
// that determines which codon it can bind to during protein synthesis. The anticodon
// forms Watson-Crick base pairs with the codon in an antiparallel fashion:
//   - Codon runs 5' → 3' on mRNA
//   - Anticodon runs 3' → 5' on tRNA (antiparallel binding)
//
// The anticodons returned here are stored in their conventional 5' → 3' reading direction
// as they appear in the tRNA structure. For example:
//   - Codon ATG (5'-ATG-3') pairs with anticodon CAU (5'-CAU-3')
//   - When aligned: 5'-ATG-3' (codon)
//                   3'-GUA-5' (anticodon, antiparallel)
//                   5'-CAU-3' (same anticodon, conventional notation)
std::string TRNA::getAnticodon(StringDict::ID tRNAId)
{
    switch (tRNAId)
    {
        // Start codon
        case StringDict::ID::TRNA_MET_ATG:
        case StringDict::ID::TRNA_MET_ATG_CHARGED:        return "CAU";
        
        // Common amino acids  
        case StringDict::ID::TRNA_GLY_GGA:
        case StringDict::ID::TRNA_GLY_GGA_CHARGED:        return "UCC";
        case StringDict::ID::TRNA_GLY_GGT:
        case StringDict::ID::TRNA_GLY_GGT_CHARGED:        return "ACC";
        case StringDict::ID::TRNA_ALA_GCA:
        case StringDict::ID::TRNA_ALA_GCA_CHARGED:        return "UGC";
        case StringDict::ID::TRNA_ALA_GCC:
        case StringDict::ID::TRNA_ALA_GCC_CHARGED:        return "GGC";
        case StringDict::ID::TRNA_LEU_CTG:
        case StringDict::ID::TRNA_LEU_CTG_CHARGED:        return "CAG";
        case StringDict::ID::TRNA_LEU_CTC:
        case StringDict::ID::TRNA_LEU_CTC_CHARGED:        return "GAG";
        case StringDict::ID::TRNA_SER_TCA:
        case StringDict::ID::TRNA_SER_TCA_CHARGED:        return "UGA";
        case StringDict::ID::TRNA_SER_TCG:
        case StringDict::ID::TRNA_SER_TCG_CHARGED:        return "CGA";
        case StringDict::ID::TRNA_VAL_GTG:
        case StringDict::ID::TRNA_VAL_GTG_CHARGED:        return "CAC";
        case StringDict::ID::TRNA_VAL_GTC:
        case StringDict::ID::TRNA_VAL_GTC_CHARGED:        return "GAC";
        
        // Less common but essential amino acids
        case StringDict::ID::TRNA_PRO_CCA:
        case StringDict::ID::TRNA_PRO_CCA_CHARGED:        return "UGG";
        case StringDict::ID::TRNA_THR_ACA:
        case StringDict::ID::TRNA_THR_ACA_CHARGED:        return "GGU";
        case StringDict::ID::TRNA_ASP_GAC:
        case StringDict::ID::TRNA_ASP_GAC_CHARGED:        return "GUC";
        case StringDict::ID::TRNA_GLU_GAG:
        case StringDict::ID::TRNA_GLU_GAG_CHARGED:        return "CUC";
        case StringDict::ID::TRNA_LYS_AAG:
        case StringDict::ID::TRNA_LYS_AAG_CHARGED:        return "CUU";
        case StringDict::ID::TRNA_ARG_CGA:
        case StringDict::ID::TRNA_ARG_CGA_CHARGED:        return "UCG";
        case StringDict::ID::TRNA_HIS_CAC:
        case StringDict::ID::TRNA_HIS_CAC_CHARGED:        return "GUG";
        case StringDict::ID::TRNA_PHE_TTC:
        case StringDict::ID::TRNA_PHE_TTC_CHARGED:        return "GAA";
        case StringDict::ID::TRNA_TYR_TAC:
        case StringDict::ID::TRNA_TYR_TAC_CHARGED:        return "GUA";
        case StringDict::ID::TRNA_CYS_TGC:
        case StringDict::ID::TRNA_CYS_TGC_CHARGED:        return "GCA";
        case StringDict::ID::TRNA_TRP_TGG:
        case StringDict::ID::TRNA_TRP_TGG_CHARGED:        return "CCA";
        case StringDict::ID::TRNA_ASN_AAC:
        case StringDict::ID::TRNA_ASN_AAC_CHARGED:        return "GUU";
        case StringDict::ID::TRNA_GLN_CAG:
        case StringDict::ID::TRNA_GLN_CAG_CHARGED:        return "CUG";
        case StringDict::ID::TRNA_ILE_ATC:
        case StringDict::ID::TRNA_ILE_ATC_CHARGED:        return "GAU";
        
        default:
            return ""; // Unknown tRNA ID
    }
}

const std::array<StringDict::ID, 25>& TRNA::getUnchargedTRNAIds()
{
    static const std::array<StringDict::ID, 25> unchargedTRNAIds = {{
        // Start codon
        StringDict::ID::TRNA_MET_ATG,
        
        // Common amino acids
        StringDict::ID::TRNA_GLY_GGA,
        StringDict::ID::TRNA_GLY_GGT,
        StringDict::ID::TRNA_ALA_GCA,
        StringDict::ID::TRNA_ALA_GCC,
        StringDict::ID::TRNA_LEU_CTG,
        StringDict::ID::TRNA_LEU_CTC,
        StringDict::ID::TRNA_SER_TCA,
        StringDict::ID::TRNA_SER_TCG,
        StringDict::ID::TRNA_VAL_GTG,
        StringDict::ID::TRNA_VAL_GTC,
        
        // Less common but essential amino acids
        StringDict::ID::TRNA_PRO_CCA,
        StringDict::ID::TRNA_THR_ACA,
        StringDict::ID::TRNA_ASP_GAC,
        StringDict::ID::TRNA_GLU_GAG,
        StringDict::ID::TRNA_LYS_AAG,
        StringDict::ID::TRNA_ARG_CGA,
        StringDict::ID::TRNA_HIS_CAC,
        StringDict::ID::TRNA_PHE_TTC,
        StringDict::ID::TRNA_TYR_TAC,
        StringDict::ID::TRNA_CYS_TGC,
        StringDict::ID::TRNA_TRP_TGG,
        StringDict::ID::TRNA_ASN_AAC,
        StringDict::ID::TRNA_GLN_CAG,
        StringDict::ID::TRNA_ILE_ATC
    }};
    return unchargedTRNAIds;
}

// Converts a codon sequence to its corresponding anticodon sequence
//
// BIOLOGICAL BACKGROUND:
// During protein synthesis, tRNA anticodons bind to mRNA codons through Watson-Crick
// base pairing in an antiparallel orientation. This function simulates that process:
//
// 1. BASE PAIRING RULES:
//    - A (adenine) pairs with U (uracil) 
//    - T (thymine, DNA only) pairs with A (adenine)
//    - U (uracil, RNA) pairs with A (adenine)
//    - G (guanine) pairs with C (cytosine)
//    - C (cytosine) pairs with G (guanine)
//
// 2. ANTIPARALLEL BINDING:
//    - Codon:     5'-ATG-3' (on mRNA, read left to right)
//    - Anticodon: 3'-UAC-5' (on tRNA, binds antiparallel)
//    - But conventionally written as: 5'-CAU-3' (reversed)
//
// 3. WHY REVERSE?
//    The anticodon binds in the opposite direction to the codon. To get the
//    conventional 5' → 3' notation of the anticodon, we must reverse the
//    complementary sequence. This ensures compatibility with stored tRNA
//    anticodon sequences in getAnticodon().
//
// EXAMPLE: ATG → complement UAC → reverse CAU (final anticodon)
std::string TRNA::codonToAnticodon(const std::string& codon)
{
    if (codon.length() != 3) return "";
    
    std::string anticodon = codon;
    for (char& c : anticodon)
    {
        switch (c)
        {
            case 'A': c = 'U'; break;  // A pairs with U
            case 'T': c = 'A'; break;  // T (DNA) pairs with A  
            case 'U': c = 'A'; break;  // U (RNA) pairs with A
            case 'G': c = 'C'; break;  // G pairs with C
            case 'C': c = 'G'; break;  // C pairs with G
            default: return ""; // Invalid nucleotide
        }
    }
    
    // Reverse the anticodon (anticodons read 3' to 5' opposite to codons)
    std::reverse(anticodon.begin(), anticodon.end());
    return anticodon;
}

std::vector<StringDict::ID> TRNA::getChargedTRNAsWithAnticodon(const std::string& anticodon)
{
    std::vector<StringDict::ID> matchingTRNAs;
    
    // Check all charged tRNA IDs to see which ones have this anticodon
    static const StringDict::ID chargedTRNAIds[] = {
        StringDict::ID::TRNA_MET_ATG_CHARGED,
        StringDict::ID::TRNA_GLY_GGA_CHARGED,
        StringDict::ID::TRNA_GLY_GGT_CHARGED,
        StringDict::ID::TRNA_ALA_GCA_CHARGED,
        StringDict::ID::TRNA_ALA_GCC_CHARGED,
        StringDict::ID::TRNA_LEU_CTG_CHARGED,
        StringDict::ID::TRNA_LEU_CTC_CHARGED,
        StringDict::ID::TRNA_SER_TCA_CHARGED,
        StringDict::ID::TRNA_SER_TCG_CHARGED,
        StringDict::ID::TRNA_VAL_GTG_CHARGED,
        StringDict::ID::TRNA_VAL_GTC_CHARGED,
        StringDict::ID::TRNA_PRO_CCA_CHARGED,
        StringDict::ID::TRNA_THR_ACA_CHARGED,
        StringDict::ID::TRNA_ASP_GAC_CHARGED,
        StringDict::ID::TRNA_GLU_GAG_CHARGED,
        StringDict::ID::TRNA_LYS_AAG_CHARGED,
        StringDict::ID::TRNA_ARG_CGA_CHARGED,
        StringDict::ID::TRNA_HIS_CAC_CHARGED,
        StringDict::ID::TRNA_PHE_TTC_CHARGED,
        StringDict::ID::TRNA_TYR_TAC_CHARGED,
        StringDict::ID::TRNA_CYS_TGC_CHARGED,
        StringDict::ID::TRNA_TRP_TGG_CHARGED,
        StringDict::ID::TRNA_ASN_AAC_CHARGED,
        StringDict::ID::TRNA_GLN_CAG_CHARGED,
        StringDict::ID::TRNA_ILE_ATC_CHARGED
    };
    
    for (StringDict::ID chargedID : chargedTRNAIds) {
        if (getAnticodon(chargedID) == anticodon) {
            matchingTRNAs.push_back(chargedID);
        }
    }
    
    return matchingTRNAs;
}

void TRNA::runTests()
{
    // Test isChargedTRNA function
    // Positive tests - loop through all charged tRNA IDs
    constexpr int EXPECTED_CHARGED_TRNA_COUNT = 25;  // Total number of charged tRNA variants
    int chargedTRNACount = 0;
    
    for (int i = static_cast<int>(StringDict::ID::TRNA_MET_ATG_CHARGED); i <= static_cast<int>(StringDict::ID::TRNA_ILE_ATC_CHARGED); ++i) {
        StringDict::ID currentID = static_cast<StringDict::ID>(i);
        chargedTRNACount++;
        
        if (!isChargedTRNA(currentID)) {
            const std::string& name = StringDict::idToString(currentID);
            LOG_ERROR("isChargedTRNA test failed: %s (ID %d) should be detected as charged but wasn't", name.c_str(), i);
            assert(false && "Charged tRNA not detected correctly");
        }
    }
    
    // Verify we tested the expected number of charged tRNAs
    if (chargedTRNACount != EXPECTED_CHARGED_TRNA_COUNT) {
        LOG_ERROR("isChargedTRNA test failed: Expected %d charged tRNAs but found %d", EXPECTED_CHARGED_TRNA_COUNT, chargedTRNACount);
        assert(chargedTRNACount == EXPECTED_CHARGED_TRNA_COUNT && "Unexpected number of charged tRNAs");
    }
    
    // Negative tests - these should return false
    if (isChargedTRNA(StringDict::ID::TRNA_MET_ATG)) {
        LOG_ERROR("isChargedTRNA test failed: TRNA_MET_ATG (uncharged) was incorrectly detected as charged");
        assert(false && "Uncharged tRNA incorrectly detected as charged");
    }
    if (isChargedTRNA(StringDict::ID::TRNA_GLY_GGA)) {
        LOG_ERROR("isChargedTRNA test failed: TRNA_GLY_GGA (uncharged) was incorrectly detected as charged");
        assert(false && "Uncharged tRNA incorrectly detected as charged");
    }
    if (isChargedTRNA(StringDict::ID::PAR_1)) {
        LOG_ERROR("isChargedTRNA test failed: PAR_1 (non-tRNA) was incorrectly detected as charged");
        assert(false && "Non-tRNA incorrectly detected as charged");
    }
    if (isChargedTRNA(StringDict::ID::ATP)) {
        LOG_ERROR("isChargedTRNA test failed: ATP (non-tRNA) was incorrectly detected as charged");
        assert(false && "Non-tRNA incorrectly detected as charged");
    }
    if (isChargedTRNA(StringDict::ID::eUNKNOWN)) {
        LOG_ERROR("isChargedTRNA test failed: eUNKNOWN was incorrectly detected as charged");
        assert(false && "eUNKNOWN incorrectly detected as charged");
    }
    
    // Test getChargedVariant function
    // Test a few representative cases
    if (getChargedVariant(StringDict::ID::TRNA_MET_ATG) != StringDict::ID::TRNA_MET_ATG_CHARGED) {
        LOG_ERROR("getChargedVariant test failed: TRNA_MET_ATG should map to TRNA_MET_ATG_CHARGED");
        assert(false && "getChargedVariant mapping incorrect");
    }
    if (getChargedVariant(StringDict::ID::TRNA_GLY_GGA) != StringDict::ID::TRNA_GLY_GGA_CHARGED) {
        LOG_ERROR("getChargedVariant test failed: TRNA_GLY_GGA should map to TRNA_GLY_GGA_CHARGED");
        assert(false && "getChargedVariant mapping incorrect");
    }
    if (getChargedVariant(StringDict::ID::TRNA_ILE_ATC) != StringDict::ID::TRNA_ILE_ATC_CHARGED) {
        LOG_ERROR("getChargedVariant test failed: TRNA_ILE_ATC should map to TRNA_ILE_ATC_CHARGED");
        assert(false && "getChargedVariant mapping incorrect");
    }
    
    // Test codonToAnticodon function
    // Test DNA codon (with T) - ATG complement is UAC, reversed is CAU
    if (codonToAnticodon("ATG") != "CAU") {
        LOG_ERROR("codonToAnticodon test failed: ATG (DNA) should convert to CAU");
        assert(false && "codonToAnticodon DNA conversion incorrect");
    }
    // Test RNA codon (with U) - AUG complement is UAC, reversed is CAU
    if (codonToAnticodon("AUG") != "CAU") {
        LOG_ERROR("codonToAnticodon test failed: AUG (RNA) should convert to CAU");
        assert(false && "codonToAnticodon RNA conversion incorrect");
    }
    // Test GGA - complement is CCU, reversed is UCC
    if (codonToAnticodon("GGA") != "UCC") {
        LOG_ERROR("codonToAnticodon test failed: GGA should convert to UCC");
        assert(false && "codonToAnticodon conversion incorrect");
    }
    
    // Test getChargedTRNAsWithAnticodon function
    auto metTRNAs = getChargedTRNAsWithAnticodon("CAU");
    if (std::find(metTRNAs.begin(), metTRNAs.end(), StringDict::ID::TRNA_MET_ATG_CHARGED) == metTRNAs.end()) {
        LOG_ERROR("getChargedTRNAsWithAnticodon test failed: CAU anticodon should match TRNA_MET_ATG_CHARGED");
        assert(false && "getChargedTRNAsWithAnticodon lookup incorrect");
    }
    
    // Test that GGA codon converts to correct anticodon and finds matching tRNA
    auto glycineTRNAs = getChargedTRNAsWithAnticodon("UCC");
    if (std::find(glycineTRNAs.begin(), glycineTRNAs.end(), StringDict::ID::TRNA_GLY_GGA_CHARGED) == glycineTRNAs.end()) {
        LOG_ERROR("getChargedTRNAsWithAnticodon test failed: UCC anticodon should match TRNA_GLY_GGA_CHARGED");
        assert(false && "getChargedTRNAsWithAnticodon lookup for GGA incorrect");
    }
}
