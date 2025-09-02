#include "pch.h"
#include "TRNA.h"
#include <algorithm>
#include <sstream>

void TRNA::charge(double dt)
{
    if (!m_bCharged)
    {
        // Probability of charging in this time step
        double chargeProbability = m_fChargingRate * dt;
        if ((double)rand() / RAND_MAX < chargeProbability)
        {
            m_bCharged = true;
        }
    }
}

void TRNA::discharge()
{
    m_bCharged = false;
}

std::string TRNA::getAnticodon() const
{
    switch (m_id)
    {
        // Start codon
        case StringDict::ID::TRNA_MET_ATG:        return "CAU";
        
        // Common amino acids  
        case StringDict::ID::TRNA_GLY_GGA:        return "CCU";
        case StringDict::ID::TRNA_GLY_GGT:        return "CCA";
        case StringDict::ID::TRNA_ALA_GCA:        return "CGU";
        case StringDict::ID::TRNA_ALA_GCC:        return "CGG";
        case StringDict::ID::TRNA_LEU_CTG:        return "GAC";
        case StringDict::ID::TRNA_LEU_CTC:        return "GAG";
        case StringDict::ID::TRNA_SER_TCA:        return "AGA";
        case StringDict::ID::TRNA_SER_TCG:        return "AGC";
        case StringDict::ID::TRNA_VAL_GTG:        return "GAC";
        case StringDict::ID::TRNA_VAL_GTC:        return "GAG";
        
        // Less common but essential amino acids
        case StringDict::ID::TRNA_PRO_CCA:        return "GGU";
        case StringDict::ID::TRNA_THR_ACA:        return "GGU";
        case StringDict::ID::TRNA_ASP_GAC:        return "CUG";
        case StringDict::ID::TRNA_GLU_GAG:        return "CUC";
        case StringDict::ID::TRNA_LYS_AAG:        return "UUC";
        case StringDict::ID::TRNA_ARG_CGA:        return "GCU";
        case StringDict::ID::TRNA_HIS_CAC:        return "GUG";
        case StringDict::ID::TRNA_PHE_TTC:        return "AAG";
        case StringDict::ID::TRNA_TYR_TAC:        return "AUG";
        case StringDict::ID::TRNA_CYS_TGC:        return "GCA";
        case StringDict::ID::TRNA_TRP_TGG:        return "CCA";
        case StringDict::ID::TRNA_ASN_AAC:        return "GUU";
        case StringDict::ID::TRNA_GLN_CAG:        return "GUC";
        case StringDict::ID::TRNA_ILE_ATC:        return "GAU";
        
        default:
            assert(false && "Unknown tRNA ID - add case to getAnticodon() switch statement");
            return ""; // Unknown tRNA ID
    }
}

bool TRNA::matchesCodon(const std::string& codon) const
{
    // Convert codon to anticodon (simplified - in reality this would be more complex)
    std::string complementary = codon;
    for (char& c : complementary)
    {
        switch (c)
        {
            case 'A': c = 'U'; break;
            case 'U': c = 'A'; break;
            case 'G': c = 'C'; break;
            case 'C': c = 'G'; break;
        }
    }
    return complementary == getAnticodon();
}
