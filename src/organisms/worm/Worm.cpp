#include "pch.h"
#include "Worm.h"
#include "chemistry/StringDict.h"
#include "biology/organelles/Cell.h"
#include "biology/simulation/CellSim.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/CellTypes.h"
#include "chemistry/Molecule.h"
#include "chemistry/TRNA.h"
#include "biology/organelles/Medium.h"
#include "biology/organelles/Cortex.h"
#include "biology/organelles/Spindle.h"
#include "chemistry/ProteinWiki.h"
#include "utils/log/ILog.h"
#include <chrono> // For high_resolution_clock
#include <cmath>  // For std::abs

// Validation thresholds based on experimental data
static constexpr double ANTERIOR_POSTERIOR_RATIO_THRESHOLD = 3.0;  // Minimum ratio for proper PAR polarization
static constexpr double NUCLEAR_SIZE_THRESHOLD = 0.8;             // Relative to initial size
static constexpr double ASYMMETRIC_DIVISION_RATIO = 0.6;         // Ratio of anterior to posterior cell size

// Development timing constants (in seconds)
static constexpr float POLARITY_ESTABLISHMENT_END_SEC = 360.0f;    // 6 minutes
static constexpr float POLARITY_MAINTENANCE_END_SEC = 600.0f;      // 10 minutes
static constexpr float NUCLEAR_ENVELOPE_BREAKDOWN_SEC = 750.0f;    // 12.5 minutes
static constexpr float SPINDLE_ASSEMBLY_START_SEC = 900.0f;        // 15 minutes
static constexpr float DIVISION_START_SEC = 1100.0f;               // 18.3 minutes

std::vector<Chromosome> Worm::initializeGenes()
{
    // Create chromosomes (C. elegans has 6 chromosomes)
    std::vector<Chromosome> chromosomes;
    chromosomes.reserve(6);

    // Create DNA for each chromosome
    auto pDNA1 = std::make_shared<DNA>();  // Chromosome I
    auto pDNA2 = std::make_shared<DNA>();  // Chromosome II
    auto pDNA3 = std::make_shared<DNA>();  // Chromosome III
    auto pDNA4 = std::make_shared<DNA>();  // Chromosome IV
    auto pDNA5 = std::make_shared<DNA>();  // Chromosome V
    auto pDNA6 = std::make_shared<DNA>();  // Chromosome X

    // Distribute genes across chromosomes (based on C. elegans genome)
    // Chromosome I - Cell fate and tRNA genes
    pDNA1->addGene(StringDict::ID::MEX_3, 0.8, 0.1);  // Anterior fate
    pDNA1->addGene(StringDict::ID::PLK_1, 1.2, 0.2);  // Polo-like kinase
    
    // Essential start codon tRNA (high abundance needed)
    pDNA1->addGene(StringDict::ID::TRNA_MET_ATG, 1.5, 0.3);  // Methionine initiator
    
    // Common amino acid tRNAs
    pDNA1->addGene(StringDict::ID::TRNA_GLY_GGA, 1.2, 0.2);  // Preferred Gly codon
    pDNA1->addGene(StringDict::ID::TRNA_ALA_GCA, 1.0, 0.15); // Preferred Ala codon
    pDNA1->addGene(StringDict::ID::TRNA_LEU_CTG, 1.4, 0.25); // Highly preferred Leu codon

    // Chromosome II - Cell cycle and tRNA genes
    pDNA2->addGene(StringDict::ID::SKN_1, 0.8, 0.1);  // Endoderm specification
    pDNA2->addGene(StringDict::ID::CYB_1, 1.2, 0.2);  // Cyclin B
    
    // More common amino acid tRNAs
    pDNA2->addGene(StringDict::ID::TRNA_GLY_GGT, 0.8, 0.15); // Second Gly choice
    pDNA2->addGene(StringDict::ID::TRNA_ALA_GCC, 0.7, 0.12); // Second Ala choice
    pDNA2->addGene(StringDict::ID::TRNA_LEU_CTC, 0.9, 0.18); // Second Leu choice
    pDNA2->addGene(StringDict::ID::TRNA_SER_TCA, 1.1, 0.2);  // Common Ser codon
    pDNA2->addGene(StringDict::ID::TRNA_VAL_GTG, 1.0, 0.18); // Preferred Val codon

    // Chromosome III - Cell cycle and tRNA genes
    pDNA3->addGene(StringDict::ID::PAL_1, 0.8, 0.1);  // Posterior fate
    pDNA3->addGene(StringDict::ID::CDK_1, 1.2, 0.2);  // Cell cycle control
    
    // Additional tRNAs for abundant amino acids
    pDNA3->addGene(StringDict::ID::TRNA_SER_TCG, 0.8, 0.15); // Alternative Ser
    pDNA3->addGene(StringDict::ID::TRNA_VAL_GTC, 0.7, 0.14); // Alternative Val
    pDNA3->addGene(StringDict::ID::TRNA_PRO_CCA, 0.9, 0.16); // Proline
    pDNA3->addGene(StringDict::ID::TRNA_THR_ACA, 0.9, 0.16); // Threonine
    pDNA3->addGene(StringDict::ID::TRNA_ASP_GAC, 1.0, 0.18); // Aspartic acid

    // Chromosome IV - Germline and charged amino acid tRNAs
    pDNA4->addGene(StringDict::ID::PIE_1, 0.8, 0.1);  // Germline specification
    
    // Charged amino acids (important for protein structure)
    pDNA4->addGene(StringDict::ID::TRNA_GLU_GAG, 1.0, 0.18); // Glutamic acid
    pDNA4->addGene(StringDict::ID::TRNA_LYS_AAG, 1.1, 0.2);  // Lysine
    pDNA4->addGene(StringDict::ID::TRNA_ARG_CGA, 0.8, 0.15); // Arginine
    pDNA4->addGene(StringDict::ID::TRNA_HIS_CAC, 0.7, 0.13); // Histidine
    pDNA4->addGene(StringDict::ID::TRNA_ASN_AAC, 0.8, 0.15); // Asparagine

    // Chromosome V - Centrosome, cytoskeleton, and aromatic amino acid tRNAs
    pDNA5->addGene(StringDict::ID::GAMMA_TUBULIN, 0.1, 0.05);  // γ-tubulin: low basal expression, will be regulated by CDK2/CyclinE
    
    // Aromatic and special amino acids (lower abundance)
    pDNA5->addGene(StringDict::ID::TRNA_PHE_TTC, 0.8, 0.14); // Phenylalanine
    pDNA5->addGene(StringDict::ID::TRNA_TYR_TAC, 0.7, 0.12); // Tyrosine
    pDNA5->addGene(StringDict::ID::TRNA_TRP_TGG, 0.5, 0.08); // Tryptophan (rare)
    pDNA5->addGene(StringDict::ID::TRNA_CYS_TGC, 0.6, 0.1);  // Cysteine
    pDNA5->addGene(StringDict::ID::TRNA_GLN_CAG, 0.9, 0.16); // Glutamine
    pDNA5->addGene(StringDict::ID::TRNA_ILE_ATC, 0.8, 0.15); // Isoleucine

    // Create chromosomes with their respective DNA
    chromosomes.emplace_back(pDNA1);
    chromosomes.emplace_back(pDNA2);
    chromosomes.emplace_back(pDNA3);
    chromosomes.emplace_back(pDNA4);
    chromosomes.emplace_back(pDNA5);  // Centrosome genes
    chromosomes.emplace_back(pDNA6);  // Empty for now

    return chromosomes;
}

void Worm::addMaternalTRNAs(Medium& medium, const float3& position)
{
    // Add essential maternal tRNAs to bootstrap translation
    // These represent tRNAs inherited from the mother egg
    
    // START CODON - absolutely essential for translation initiation
    auto pMet = std::make_shared<TRNA>("Met", "CAU", 500.0, 0.8);  // High abundance for start codon
    pMet->charge(1.0);  // Start fully charged
    medium.addTRNA(pMet, position);
    
    // MOST COMMON AMINO ACIDS - needed for early protein synthesis
    // Glycine (highly abundant in C. elegans)
    auto pGlyGGA = std::make_shared<TRNA>("Gly", "CCU", 300.0, 0.9);
    pGlyGGA->charge(1.0);
    medium.addTRNA(pGlyGGA, position);
    
    auto pGlyGGT = std::make_shared<TRNA>("Gly", "CCA", 200.0, 0.8);
    pGlyGGT->charge(1.0);
    medium.addTRNA(pGlyGGT, position);
    
    // Alanine (very common)
    auto pAlaGCA = std::make_shared<TRNA>("Ala", "CGU", 250.0, 0.9);
    pAlaGCA->charge(1.0);
    medium.addTRNA(pAlaGCA, position);
    
    auto pAlaGCC = std::make_shared<TRNA>("Ala", "CGG", 150.0, 0.8);
    pAlaGCC->charge(1.0);
    medium.addTRNA(pAlaGCC, position);
    
    // Leucine (highly preferred in C. elegans)
    auto pLeuCTG = std::make_shared<TRNA>("Leu", "GAC", 350.0, 0.9);
    pLeuCTG->charge(1.0);
    medium.addTRNA(pLeuCTG, position);
    
    auto pLeuCTC = std::make_shared<TRNA>("Leu", "GAG", 200.0, 0.8);
    pLeuCTC->charge(1.0);
    medium.addTRNA(pLeuCTC, position);
    
    // Serine (common)
    auto pSerTCA = std::make_shared<TRNA>("Ser", "AGA", 220.0, 0.8);
    pSerTCA->charge(1.0);
    medium.addTRNA(pSerTCA, position);
    
    auto pSerTCG = std::make_shared<TRNA>("Ser", "AGC", 150.0, 0.7);
    pSerTCG->charge(1.0);
    medium.addTRNA(pSerTCG, position);
    
    // Valine
    auto pValGTG = std::make_shared<TRNA>("Val", "GAC", 200.0, 0.8);
    pValGTG->charge(1.0);
    medium.addTRNA(pValGTG, position);
    
    // ESSENTIAL AMINO ACIDS - lower abundance but necessary
    // Lysine (positively charged, important for proteins)
    auto pLysAAG = std::make_shared<TRNA>("Lys", "UUC", 180.0, 0.8);
    pLysAAG->charge(1.0);
    medium.addTRNA(pLysAAG, position);
    
    // Aspartic acid (negatively charged)
    auto pAspGAC = std::make_shared<TRNA>("Asp", "CUG", 160.0, 0.8);
    pAspGAC->charge(1.0);
    medium.addTRNA(pAspGAC, position);
    
    // Glutamic acid (negatively charged)
    auto pGluGAG = std::make_shared<TRNA>("Glu", "CUC", 170.0, 0.8);
    pGluGAG->charge(1.0);
    medium.addTRNA(pGluGAG, position);
    
    // Proline (structure-forming)
    auto pProCCA = std::make_shared<TRNA>("Pro", "GGU", 140.0, 0.7);
    pProCCA->charge(1.0);
    medium.addTRNA(pProCCA, position);
    
    // Threonine
    auto pThrACA = std::make_shared<TRNA>("Thr", "GGU", 140.0, 0.7);
    pThrACA->charge(1.0);
    medium.addTRNA(pThrACA, position);
    
    // These maternal tRNAs will allow initial translation of:
    // 1. More tRNAs (from transcribed tRNA mRNAs)
    // 2. tRNA synthetases (to charge more tRNAs)
    // 3. Ribosomal proteins (to make more ribosomes)
    // 4. Other essential proteins for cellular function
}

std::shared_ptr<Medium> Worm::createZygoteMedium()
{
    // Create the internal medium
    std::shared_ptr<Medium> pInternalMedium = std::make_shared<Medium>();

    // Create and add anterior proteins at the anterior cortex
    MPopulation par3(StringDict::idToString(StringDict::ID::PAR_3), 3.9e5);
    pInternalMedium->addProtein(par3, float3(0, 1.f, 0));

    MPopulation par6(StringDict::idToString(StringDict::ID::PAR_6), 3.9e5);
    pInternalMedium->addProtein(par6, float3(0, 1.f, 0));

    MPopulation pkc3(StringDict::idToString(StringDict::ID::PKC_3), 3.9e5);
    pInternalMedium->addProtein(pkc3, float3(0, 1.f, 0));

    // Create and add posterior proteins at the posterior cortex
    MPopulation par1(StringDict::idToString(StringDict::ID::PAR_1), 3.9e5);
    pInternalMedium->addProtein(par1, float3(0, -1.f, 0));

    MPopulation par2(StringDict::idToString(StringDict::ID::PAR_2), 3.9e5);
    pInternalMedium->addProtein(par2, float3(0, -1.f, 0));

    // Initialize maternal proteins at cell center
    float3 center(0.0f, 0.0f, 0.0f);

    // Add maternal CDK-1 and CYB-1 (Cyclin B)
    MPopulation cdk1(StringDict::idToString(StringDict::ID::CDK_1), 1500.0);  // Initial amount above threshold (1000)
    MPopulation cyb1(StringDict::idToString(StringDict::ID::CYB_1), 1500.0);  // Initial amount above threshold (1000)
    pInternalMedium->addProtein(cdk1, center);
    pInternalMedium->addProtein(cyb1, center);

    // Add centrosome-related proteins for proper centrosome function
    MPopulation cdk2(StringDict::idToString(StringDict::ID::CDK_2), 800.0);  // For centrosome duplication
    MPopulation cce1(StringDict::idToString(StringDict::ID::CCE_1), 800.0);  // Cyclin E for centrosome duplication
    pInternalMedium->addProtein(cdk2, center);
    pInternalMedium->addProtein(cce1, center);
    
    // Add maternal ATP for translation
    pInternalMedium->addATP(50000.0, center);  // Sufficient ATP for early translation
    
    // Add maternal tRNAs (essential for translation bootstrap)
    // Without these, mRNAs (including tRNA mRNAs) cannot be translated
    addMaternalTRNAs(*pInternalMedium, center);

    return pInternalMedium;
}

Worm::Worm()
{
    // Initialize the string dictionary first
    StringDict::initialize();
    
    auto chromosomes = initializeGenes();
    auto pInternalMedium = createZygoteMedium();
    auto pCell = Cell::createCell(pInternalMedium, chromosomes);
    
    // Simulate fertilization by adding a centrosome from the sperm
    // In C. elegans, the egg lacks a centrosome and receives one from the sperm
    // The sperm enters at the posterior end, so place centrosome near posterior boundary
    float3 posteriorEntryPoint(0.0f, -0.8f, 0.0f);  // Near posterior boundary
    auto pCentrosome = std::make_shared<Centrosome>(std::weak_ptr<Cell>(pCell), posteriorEntryPoint);
    pCell->addOrganelle(StringDict::ID::ORGANELLE_CENTROSOME, pCentrosome);
    
    // Create CellSim to wrap the Cell
    auto pCellSim = std::make_shared<CellSim>(pCell);
    m_pCellSims.push_back(pCellSim);
    
    // Set up the data collector
    setupDataCollector();
}

void Worm::setupDataCollector()
{
    // Get access to the cell's internal medium
    if (m_pCellSims.empty()) {
        LOG_ERROR("Cannot set up data collector: no cells available");
        return;
    }
    
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    
    // Initialize the DataCollector with the cell's internal medium
    m_pDataCollector = std::make_unique<DataCollector>(
        internalMedium, 
        "worm_simulation_data.csv",
        0.1f  // Collect data every 0.1 seconds
    );
    
    // Add collection points for anterior and posterior positions
    float3 anteriorPos(0.0f, 1.f, 0.0f);  // Anterior position
    float3 posteriorPos(0.0f, -1.f, 0.0f); // Posterior position
    float3 centerPos(0.0f, 0.0f, 0.0f);   // Center position for centrosome tracking
    
    // Get membrane-bound protein names using the utility function
    std::string par2Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::ORGANELLE_CORTEX);
    std::string par3Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::ORGANELLE_CORTEX);
    std::string bindingSiteCortex = StringDict::idToString(StringDict::ID::ORGANELLE_CORTEX);

    // Add collection points with specific proteins to track
    m_pDataCollector->addCollectionPoint(
        posteriorPos,
        "Posterior", 
        { StringDict::idToString(StringDict::ID::GAMMA_TUBULIN) }
    );
}

void Worm::simulateStep(const TimeContext& time)
{
    // Use local variables for timing
    auto stepStartTime = std::chrono::high_resolution_clock::now();
    
    // Call the base class simulateStep to handle standard simulation
    Organism::simulateStep(time);
    
    // Calculate time taken for this step using local variables
    auto stepEndTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stepEndTime - stepStartTime);
    double stepTimeMsec = duration.count() / 1000.0; // Convert to milliseconds
    
    // Update total simulation time
    m_fTotalTime += static_cast<float>(time.m_deltaTSec);
    
    // Force data collection at the current time (regardless of interval)
    if (m_pDataCollector) {
        m_pDataCollector->forceCollection(m_fTotalTime, stepTimeMsec);
    }
}

bool Worm::validatePARPolarization(float fTimeSec) const
{
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    
    float3 anteriorPos(0.0f, 1.f, 0.0f);
    float3 posteriorPos(0.0f, -1.f, 0.0f);
    
    // Get membrane-bound protein names using the utility function
    std::string par3Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::ORGANELLE_CORTEX);
    std::string par2Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::ORGANELLE_CORTEX);
    
    // Check membrane-bound proteins
    double anteriorPAR3 = internalMedium.getProteinNumber(par3Membrane, anteriorPos);
    double posteriorPAR3 = internalMedium.getProteinNumber(par3Membrane, posteriorPos);
    double anteriorPAR2 = internalMedium.getProteinNumber(par2Membrane, anteriorPos);
    double posteriorPAR2 = internalMedium.getProteinNumber(par2Membrane, posteriorPos);
    
    // Check during polarity establishment (0-6 minutes)
    if (fTimeSec < POLARITY_ESTABLISHMENT_END_SEC) {
        if (anteriorPAR3 / (posteriorPAR3 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient anterior %s polarization at %.2lf sec", par3Membrane.c_str(), fTimeSec);
            return false;
        }
        
        if (posteriorPAR2 / (anteriorPAR2 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient posterior %s polarization at %.2lf sec", par2Membrane.c_str(), fTimeSec);
            return false;
        }
    }
    
    return true;
}

bool Worm::validateCellCycle(float fTimeSec) const
{
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    float3 nuclearPos(0.0f, 0.0f, 0.0f);
    double cdk1Level = internalMedium.getProteinNumber(StringDict::idToString(StringDict::ID::CDK_1), nuclearPos);

    // Before nuclear envelope breakdown (0-12.5 minutes): CDK-1 should be relatively low
    if (fTimeSec < NUCLEAR_ENVELOPE_BREAKDOWN_SEC && cdk1Level > 1000) {
        LOG_INFO("Warning: CDK-1 levels too high before NEBD at %.2lf sec", fTimeSec);
        return false;
    }
    
    // During mitotic entry (12.5-15 minutes): CDK-1 should increase
    if (fTimeSec >= NUCLEAR_ENVELOPE_BREAKDOWN_SEC && fTimeSec < SPINDLE_ASSEMBLY_START_SEC && cdk1Level < 1500) {
        LOG_INFO("Warning: CDK-1 levels too low during mitotic entry at %.2lf sec", fTimeSec);
        return false;
    }
    
    return true;
}

bool Worm::validateAsymmetricDivision(float fTimeSec) const
{
    // Only check during late stages (after 15 minutes)
    if (fTimeSec < SPINDLE_ASSEMBLY_START_SEC) return true;
    
    auto pSpindle = std::dynamic_pointer_cast<Spindle>(m_pCellSims[0]->getCell()->getOrganelle(StringDict::ID::ORGANELLE_SPINDLE));
    float3 spindlePos = pSpindle->getPosition();
    
    if (spindlePos.y > -0.1f) {
        LOG_INFO("Warning: Spindle not properly positioned toward posterior at %.2lf sec", fTimeSec);
        return false;
    }
    
    return true;
}

bool Worm::validateCentrosomeBehavior(float fTimeSec) const
{
    auto pCentrosome = std::dynamic_pointer_cast<Centrosome>(m_pCellSims[0]->getCell()->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
    if (!pCentrosome) {
        // Before fertilization, there should be no centrosome
        if (fTimeSec < 1.0f) {  // Assume fertilization happens within 1 second
            return true;  // This is expected before fertilization
        }
        LOG_INFO("Warning: No centrosome found in cell at %.2lf sec (after expected fertilization time)", fTimeSec);
        return false;
    }
    
    // Check centrosome position during different phases
    float3 centrosomePos = pCentrosome->getNormalizedPosition();
    CellCycleState cellCycleState = m_pCellSims[0]->getCell()->getCellCycleState();
    
    switch (cellCycleState) {
        case CellCycleState::INTERPHASE:
            // Centrosome should be near the nucleus (center)
            if (std::abs(centrosomePos.x) > 0.2f || std::abs(centrosomePos.y) > 0.2f || std::abs(centrosomePos.z) > 0.2f) {
                LOG_INFO("Warning: Centrosome too far from nucleus during interphase at %.2lf sec", fTimeSec);
                return false;
            }
            break;
            
        case CellCycleState::PROPHASE:
        case CellCycleState::METAPHASE:
            // Centrosome should be duplicated and moving to poles
            if (!pCentrosome->isDuplicated()) {
                LOG_INFO("Warning: Centrosome not duplicated during mitosis at %.2lf sec", fTimeSec);
                return false;
            }
            // Check if centrosome is moving toward poles
            if (std::abs(centrosomePos.y) < 0.5f) {
                LOG_INFO("Warning: Centrosome not properly positioned at poles during mitosis at %.2lf sec", fTimeSec);
                return false;
            }
            break;
            
        case CellCycleState::ANAPHASE:
        case CellCycleState::TELOPHASE:
            // Centrosome should be at poles
            if (std::abs(centrosomePos.y) < 0.7f) {
                LOG_INFO("Warning: Centrosome not at poles during anaphase/telophase at %.2lf sec", fTimeSec);
                return false;
            }
            break;
            
        case CellCycleState::CYTOKINESIS:
            // Centrosome should be resetting for next cycle
            if (pCentrosome->isDuplicated()) {
                LOG_INFO("Warning: Centrosome still duplicated during cytokinesis at %.2lf sec", fTimeSec);
                return false;
            }
            break;
    }
    
    return true;
}
