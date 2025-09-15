#include "MoleculeWiki.h"
#include "chemistry/interactions/PhosphorylationInteraction.h"
#include "chemistry/interactions/ComplexFormationInteraction.h"
#include "chemistry/interactions/DephosphorylationInteraction.h"
#include "BindingSurface.h"
// Interactions are now managed in InteractionsWiki (separate project)
#include "utils/log/ILog.h"
#include "utils/fileUtils/fileUtils.h"
#include <algorithm>
#include <iterator>
#include <filesystem>

// Initialize static members
std::unordered_map<Molecule, MolInfo> MoleculeWiki::m_moleculesInfo;

void MoleculeWiki::Initialize()
{
    // Initialize molecule info and interactions
    m_moleculesInfo.clear();
    
    // Initialize tRNA molecule information
    initializeTRNAInfo();
    
    // Initialize mRNA molecule information
    initializeMRNAInfo();
    
    // Interactions are initialized from InteractionsWiki
}

// Interactions accessors were removed; use InteractionsWiki directly

std::string MoleculeWiki::GetBoundProteinName(const std::string& proteinName, StringDict::ID surface)
{
    return proteinName + ":" + StringDict::idToString(surface);
}

const MolInfo& MoleculeWiki::getInfo(const Molecule& molecule)
{
    auto it = m_moleculesInfo.find(molecule);
    if (it != m_moleculesInfo.end()) {
        return it->second;
    }
    
    // If not found, log error and assert
    std::string errorMsg = "Molecule info not found for: " + molecule.getName() + " (type: " + 
                          std::to_string(static_cast<int>(molecule.getType())) + ")";
    LOG_ERROR(errorMsg.c_str());
    assert(false && "Molecule information not found in MoleculeWiki - all used molecules must be initialized");
    
    // This should never be reached, but keep for compilation
    static const MolInfo defaultInfo("No information available", "", 0.0, "", 0.0, 0.0);
    return defaultInfo;
}

void MoleculeWiki::initializeTRNAInfo()
{
    // Initialize uncharged tRNA molecules with their charging rates
    
    // Start codon - absolutely essential
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_MET_ATG, ChemicalType::TRNA)] = 
        MolInfo("Methionine tRNA (uncharged)", "tRNA", 25000.0, "Start codon tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_MET_ATG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Methionine tRNA (charged)", "tRNA", 25000.0, "Start codon tRNA", 0.0, 0.0, 0.0);
    
    // Common amino acids - high charging rates
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGA, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGA (uncharged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGA (charged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGT, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGT (uncharged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGT_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glycine tRNA GGT (charged)", "tRNA", 25000.0, "Glycine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCA, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCA (uncharged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCA (charged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCC, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCC (uncharged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Alanine tRNA GCC (charged)", "tRNA", 25000.0, "Alanine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTG, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTG (uncharged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.9);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTG (charged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTC, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTC (uncharged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Leucine tRNA CTC (charged)", "tRNA", 25000.0, "Leucine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCA, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCA (uncharged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCA (charged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCG, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCG (uncharged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Serine tRNA TCG (charged)", "tRNA", 25000.0, "Serine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTG, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTG (uncharged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTG (charged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTC, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTC (uncharged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Valine tRNA GTC (charged)", "tRNA", 25000.0, "Valine tRNA", 0.0, 0.0, 0.0);
        
    // Essential amino acids - lower charging rates
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LYS_AAG, ChemicalType::TRNA)] = 
        MolInfo("Lysine tRNA AAG (uncharged)", "tRNA", 25000.0, "Lysine tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LYS_AAG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Lysine tRNA AAG (charged)", "tRNA", 25000.0, "Lysine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASP_GAC, ChemicalType::TRNA)] = 
        MolInfo("Aspartic acid tRNA GAC (uncharged)", "tRNA", 25000.0, "Aspartic acid tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASP_GAC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Aspartic acid tRNA GAC (charged)", "tRNA", 25000.0, "Aspartic acid tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLU_GAG, ChemicalType::TRNA)] = 
        MolInfo("Glutamic acid tRNA GAG (uncharged)", "tRNA", 25000.0, "Glutamic acid tRNA", 0.0, 0.0, 0.8);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLU_GAG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glutamic acid tRNA GAG (charged)", "tRNA", 25000.0, "Glutamic acid tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PRO_CCA, ChemicalType::TRNA)] = 
        MolInfo("Proline tRNA CCA (uncharged)", "tRNA", 25000.0, "Proline tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PRO_CCA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Proline tRNA CCA (charged)", "tRNA", 25000.0, "Proline tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_THR_ACA, ChemicalType::TRNA)] = 
        MolInfo("Threonine tRNA ACA (uncharged)", "tRNA", 25000.0, "Threonine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_THR_ACA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Threonine tRNA ACA (charged)", "tRNA", 25000.0, "Threonine tRNA", 0.0, 0.0, 0.0);
        
    // Add remaining tRNAs with standard rates...
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ARG_CGA, ChemicalType::TRNA)] = 
        MolInfo("Arginine tRNA CGA (uncharged)", "tRNA", 25000.0, "Arginine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ARG_CGA_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Arginine tRNA CGA (charged)", "tRNA", 25000.0, "Arginine tRNA", 0.0, 0.0, 0.0);
    
    // Missing tRNAs that were causing assertion errors
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_HIS_CAC, ChemicalType::TRNA)] = 
        MolInfo("Histidine tRNA CAC (uncharged)", "tRNA", 25000.0, "Histidine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_HIS_CAC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Histidine tRNA CAC (charged)", "tRNA", 25000.0, "Histidine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PHE_TTC, ChemicalType::TRNA)] = 
        MolInfo("Phenylalanine tRNA TTC (uncharged)", "tRNA", 25000.0, "Phenylalanine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PHE_TTC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Phenylalanine tRNA TTC (charged)", "tRNA", 25000.0, "Phenylalanine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TYR_TAC, ChemicalType::TRNA)] = 
        MolInfo("Tyrosine tRNA TAC (uncharged)", "tRNA", 25000.0, "Tyrosine tRNA", 0.0, 0.0, 0.6);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TYR_TAC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Tyrosine tRNA TAC (charged)", "tRNA", 25000.0, "Tyrosine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_CYS_TGC, ChemicalType::TRNA)] = 
        MolInfo("Cysteine tRNA TGC (uncharged)", "tRNA", 25000.0, "Cysteine tRNA", 0.0, 0.0, 0.6);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_CYS_TGC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Cysteine tRNA TGC (charged)", "tRNA", 25000.0, "Cysteine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TRP_TGG, ChemicalType::TRNA)] = 
        MolInfo("Tryptophan tRNA TGG (uncharged)", "tRNA", 25000.0, "Tryptophan tRNA", 0.0, 0.0, 0.6);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TRP_TGG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Tryptophan tRNA TGG (charged)", "tRNA", 25000.0, "Tryptophan tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASN_AAC, ChemicalType::TRNA)] = 
        MolInfo("Asparagine tRNA AAC (uncharged)", "tRNA", 25000.0, "Asparagine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASN_AAC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Asparagine tRNA AAC (charged)", "tRNA", 25000.0, "Asparagine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLN_CAG, ChemicalType::TRNA)] = 
        MolInfo("Glutamine tRNA CAG (uncharged)", "tRNA", 25000.0, "Glutamine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLN_CAG_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Glutamine tRNA CAG (charged)", "tRNA", 25000.0, "Glutamine tRNA", 0.0, 0.0, 0.0);
    
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ILE_ATC, ChemicalType::TRNA)] = 
        MolInfo("Isoleucine tRNA ATC (uncharged)", "tRNA", 25000.0, "Isoleucine tRNA", 0.0, 0.0, 0.7);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ILE_ATC_CHARGED, ChemicalType::TRNA)] = 
        MolInfo("Isoleucine tRNA ATC (charged)", "tRNA", 25000.0, "Isoleucine tRNA", 0.0, 0.0, 0.0);
}

void MoleculeWiki::initializeMRNAInfo()
{
    // Initialize mRNA molecules with their translation rates
    // Translation rates are in proteins/second/mRNA - biological range is 0.1-10.0 proteins per mRNA per second
    
    // Cell fate specification mRNAs - critical for early development, high translation rates
    m_moleculesInfo[Molecule(StringDict::ID::PIE_1, ChemicalType::MRNA)] = 
        MolInfo("PIE-1 mRNA", "mRNA", 50000.0, "Germline specification mRNA", 1800.0, 2.0, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::PAL_1, ChemicalType::MRNA)] = 
        MolInfo("PAL-1 mRNA", "mRNA", 45000.0, "Posterior fate specification mRNA", 1800.0, 1.8, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::SKN_1, ChemicalType::MRNA)] = 
        MolInfo("SKN-1 mRNA", "mRNA", 48000.0, "Endoderm specification mRNA", 1800.0, 2.2, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::MEX_3, ChemicalType::MRNA)] = 
        MolInfo("MEX-3 mRNA", "mRNA", 46000.0, "Anterior fate specification mRNA", 1800.0, 1.5, 0.0);
    
    // Cell cycle mRNAs - essential for division timing
    m_moleculesInfo[Molecule(StringDict::ID::CDK_1, ChemicalType::MRNA)] = 
        MolInfo("CDK-1 mRNA", "mRNA", 40000.0, "Cyclin-dependent kinase mRNA", 2400.0, 3.0, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::CDK_2, ChemicalType::MRNA)] = 
        MolInfo("CDK-2 mRNA", "mRNA", 38000.0, "CDK-2 transcriptional regulator mRNA", 2400.0, 2.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::CYB_1, ChemicalType::MRNA)] = 
        MolInfo("CYB-1 mRNA", "mRNA", 42000.0, "Cyclin B mRNA", 1500.0, 2.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::CCE_1, ChemicalType::MRNA)] = 
        MolInfo("CCE-1 mRNA", "mRNA", 40000.0, "Cyclin E transcriptional regulator mRNA", 1800.0, 2.8, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::PLK_1, ChemicalType::MRNA)] = 
        MolInfo("PLK-1 mRNA", "mRNA", 38000.0, "Polo-like kinase mRNA", 2000.0, 2.8, 0.0);
    
    // Centrosome protein mRNAs - structural proteins, moderate translation
    m_moleculesInfo[Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::MRNA)] = 
        MolInfo("γ-TUBULIN mRNA", "mRNA", 44000.0, "γ-tubulin mRNA", 3600.0, 1.2, 0.0);
    
    // tRNA mRNAs - these encode the tRNA molecules themselves, low-moderate translation
    // Start codon tRNA
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_MET_ATG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Met-ATG mRNA", "mRNA", 25000.0, "Methionine tRNA gene mRNA", 7200.0, 0.8, 0.0);
    
    // Common amino acids tRNAs
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Gly-GGA mRNA", "mRNA", 25000.0, "Glycine tRNA gene mRNA", 7200.0, 0.6, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLY_GGT, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Gly-GGT mRNA", "mRNA", 25000.0, "Glycine tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Ala-GCA mRNA", "mRNA", 25000.0, "Alanine tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ALA_GCC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Ala-GCC mRNA", "mRNA", 25000.0, "Alanine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Leu-CTG mRNA", "mRNA", 25000.0, "Leucine tRNA gene mRNA", 7200.0, 0.7, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LEU_CTC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Leu-CTC mRNA", "mRNA", 25000.0, "Leucine tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Ser-TCA mRNA", "mRNA", 25000.0, "Serine tRNA gene mRNA", 7200.0, 0.6, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_SER_TCG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Ser-TCG mRNA", "mRNA", 25000.0, "Serine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Val-GTG mRNA", "mRNA", 25000.0, "Valine tRNA gene mRNA", 7200.0, 0.6, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_VAL_GTC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Val-GTC mRNA", "mRNA", 25000.0, "Valine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    
    // Less common but essential amino acids tRNAs
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PRO_CCA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Pro-CCA mRNA", "mRNA", 25000.0, "Proline tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_THR_ACA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Thr-ACA mRNA", "mRNA", 25000.0, "Threonine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASP_GAC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Asp-GAC mRNA", "mRNA", 25000.0, "Aspartic acid tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLU_GAG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Glu-GAG mRNA", "mRNA", 25000.0, "Glutamic acid tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_LYS_AAG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Lys-AAG mRNA", "mRNA", 25000.0, "Lysine tRNA gene mRNA", 7200.0, 0.5, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ARG_CGA, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Arg-CGA mRNA", "mRNA", 25000.0, "Arginine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_HIS_CAC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-His-CAC mRNA", "mRNA", 25000.0, "Histidine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_PHE_TTC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Phe-TTC mRNA", "mRNA", 25000.0, "Phenylalanine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TYR_TAC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Tyr-TAC mRNA", "mRNA", 25000.0, "Tyrosine tRNA gene mRNA", 7200.0, 0.3, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_CYS_TGC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Cys-TGC mRNA", "mRNA", 25000.0, "Cysteine tRNA gene mRNA", 7200.0, 0.3, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_TRP_TGG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Trp-TGG mRNA", "mRNA", 25000.0, "Tryptophan tRNA gene mRNA", 7200.0, 0.3, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ASN_AAC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Asn-AAC mRNA", "mRNA", 25000.0, "Asparagine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_GLN_CAG, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Gln-CAG mRNA", "mRNA", 25000.0, "Glutamine tRNA gene mRNA", 7200.0, 0.4, 0.0);
    m_moleculesInfo[Molecule(StringDict::ID::TRNA_ILE_ATC, ChemicalType::MRNA)] = 
        MolInfo("tRNA-Ile-ATC mRNA", "mRNA", 25000.0, "Isoleucine tRNA gene mRNA", 7200.0, 0.4, 0.0);
}