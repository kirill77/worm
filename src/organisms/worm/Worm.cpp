#include "pch.h"
#include "Worm.h"
#include "chemistry/molecules/StringDict.h"
#include "biology/organelles/Cell.h"
#include "biology/simulation/CellSim.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/CellTypes.h"
#include "chemistry/molecules/Molecule.h"
#include "biology/organelles/Medium.h"
#include "biology/organelles/Cortex.h"
#include "biology/organelles/Spindle.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "utils/log/ILog.h"
#include "utils/fileUtils/fileUtils.h"
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
    auto pDNA1 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome I
    auto pDNA2 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome II
    auto pDNA3 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome III
    auto pDNA4 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome IV
    auto pDNA5 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome V
    auto pDNA6 = std::make_shared<DNA>(Species::C_ELEGANS);  // Chromosome X

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
    pDNA3->addGene(StringDict::ID::CDK_2, 1.0, 0.15); // Transcriptional regulator for γ-tubulin
    pDNA3->addGene(StringDict::ID::CCE_1, 1.1, 0.18); // Cyclin E transcriptional regulator
    
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
    // We add charged tRNAs since they come from the mother in a charged state
    
    // START CODON - absolutely essential for translation initiation  
    MPopulation metCharged(Molecule(StringDict::ID::TRNA_MET_ATG_CHARGED, ChemicalType::TRNA), 500.0);
    medium.addMolecule(metCharged, position);
    
    // MOST COMMON AMINO ACIDS - needed for early protein synthesis
    // Glycine (highly abundant in C. elegans)
    MPopulation glyGGACharged(Molecule(StringDict::ID::TRNA_GLY_GGA_CHARGED, ChemicalType::TRNA), 300.0);
    medium.addMolecule(glyGGACharged, position);
    
    MPopulation glyGGTCharged(Molecule(StringDict::ID::TRNA_GLY_GGT_CHARGED, ChemicalType::TRNA), 200.0);
    medium.addMolecule(glyGGTCharged, position);
    
    // Alanine (very common)
    MPopulation alaGCACharged(Molecule(StringDict::ID::TRNA_ALA_GCA_CHARGED, ChemicalType::TRNA), 250.0);
    medium.addMolecule(alaGCACharged, position);
    
    MPopulation alaGCCCharged(Molecule(StringDict::ID::TRNA_ALA_GCC_CHARGED, ChemicalType::TRNA), 150.0);
    medium.addMolecule(alaGCCCharged, position);
    
    // Leucine (highly preferred in C. elegans)
    MPopulation leuCTGCharged(Molecule(StringDict::ID::TRNA_LEU_CTG_CHARGED, ChemicalType::TRNA), 350.0);
    medium.addMolecule(leuCTGCharged, position);
    
    MPopulation leuCTCCharged(Molecule(StringDict::ID::TRNA_LEU_CTC_CHARGED, ChemicalType::TRNA), 200.0);
    medium.addMolecule(leuCTCCharged, position);
    
    // Serine (common)
    MPopulation serTCACharged(Molecule(StringDict::ID::TRNA_SER_TCA_CHARGED, ChemicalType::TRNA), 220.0);
    medium.addMolecule(serTCACharged, position);
    
    MPopulation serTCGCharged(Molecule(StringDict::ID::TRNA_SER_TCG_CHARGED, ChemicalType::TRNA), 150.0);
    medium.addMolecule(serTCGCharged, position);
    
    // Valine (both codons needed for GAMMA_TUBULIN translation)
    MPopulation valGTGCharged(Molecule(StringDict::ID::TRNA_VAL_GTG_CHARGED, ChemicalType::TRNA), 200.0);
    medium.addMolecule(valGTGCharged, position);
    
    MPopulation valGTCCharged(Molecule(StringDict::ID::TRNA_VAL_GTC_CHARGED, ChemicalType::TRNA), 150.0);
    medium.addMolecule(valGTCCharged, position);
    
    // ESSENTIAL AMINO ACIDS - lower abundance but necessary
    // Lysine (positively charged, important for proteins)
    MPopulation lysAAGCharged(Molecule(StringDict::ID::TRNA_LYS_AAG_CHARGED, ChemicalType::TRNA), 180.0);
    medium.addMolecule(lysAAGCharged, position);
    
    // Arginine (required by γ-tubulin sequence in worm data)
    MPopulation argCGACharged(Molecule(StringDict::ID::TRNA_ARG_CGA_CHARGED, ChemicalType::TRNA), 160.0);
    medium.addMolecule(argCGACharged, position);

    // Additional essentials observed in sequence requirements
    // Histidine
    MPopulation hisCACCharged(Molecule(StringDict::ID::TRNA_HIS_CAC_CHARGED, ChemicalType::TRNA), 140.0);
    medium.addMolecule(hisCACCharged, position);
    // Tyrosine
    MPopulation tyrTACCharged(Molecule(StringDict::ID::TRNA_TYR_TAC_CHARGED, ChemicalType::TRNA), 120.0);
    medium.addMolecule(tyrTACCharged, position);
    // Cysteine
    MPopulation cysTGCCharged(Molecule(StringDict::ID::TRNA_CYS_TGC_CHARGED, ChemicalType::TRNA), 120.0);
    medium.addMolecule(cysTGCCharged, position);
    // Tryptophan
    MPopulation trpTGGCharged(Molecule(StringDict::ID::TRNA_TRP_TGG_CHARGED, ChemicalType::TRNA), 100.0);
    medium.addMolecule(trpTGGCharged, position);
    // Asparagine
    MPopulation asnAACCharged(Molecule(StringDict::ID::TRNA_ASN_AAC_CHARGED, ChemicalType::TRNA), 140.0);
    medium.addMolecule(asnAACCharged, position);
    // Glutamine
    MPopulation glnCAGCharged(Molecule(StringDict::ID::TRNA_GLN_CAG_CHARGED, ChemicalType::TRNA), 140.0);
    medium.addMolecule(glnCAGCharged, position);
    // Isoleucine
    MPopulation ileATCCharged(Molecule(StringDict::ID::TRNA_ILE_ATC_CHARGED, ChemicalType::TRNA), 150.0);
    medium.addMolecule(ileATCCharged, position);
    
    // Aspartic acid (negatively charged)
    MPopulation aspGACCharged(Molecule(StringDict::ID::TRNA_ASP_GAC_CHARGED, ChemicalType::TRNA), 160.0);
    medium.addMolecule(aspGACCharged, position);
    
    // Glutamic acid (negatively charged)
    MPopulation gluGAGCharged(Molecule(StringDict::ID::TRNA_GLU_GAG_CHARGED, ChemicalType::TRNA), 170.0);
    medium.addMolecule(gluGAGCharged, position);
    
    // Proline (structure-forming)
    MPopulation proCCACharged(Molecule(StringDict::ID::TRNA_PRO_CCA_CHARGED, ChemicalType::TRNA), 140.0);
    medium.addMolecule(proCCACharged, position);
    
    // Threonine
    MPopulation thrACACharged(Molecule(StringDict::ID::TRNA_THR_ACA_CHARGED, ChemicalType::TRNA), 140.0);
    medium.addMolecule(thrACACharged, position);
    
    // Phenylalanine (essential for GAMMA_TUBULIN translation)
    MPopulation pheTTCCharged(Molecule(StringDict::ID::TRNA_PHE_TTC_CHARGED, ChemicalType::TRNA), 130.0);
    medium.addMolecule(pheTTCCharged, position);
    
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
    MPopulation par3(Molecule(StringDict::ID::PAR_3, ChemicalType::PROTEIN, Species::C_ELEGANS), 3.9e5);
    pInternalMedium->addMolecule(par3, float3(0, 1.f, 0));

    MPopulation par6(Molecule(StringDict::ID::PAR_6, ChemicalType::PROTEIN, Species::C_ELEGANS), 3.9e5);
    pInternalMedium->addMolecule(par6, float3(0, 1.f, 0));

    MPopulation pkc3(Molecule(StringDict::ID::PKC_3, ChemicalType::PROTEIN, Species::C_ELEGANS), 3.9e5);
    pInternalMedium->addMolecule(pkc3, float3(0, 1.f, 0));

    // Create and add posterior proteins at the posterior cortex
    MPopulation par1(Molecule(StringDict::ID::PAR_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 3.9e5);
    pInternalMedium->addMolecule(par1, float3(0, -1.f, 0));

    MPopulation par2(Molecule(StringDict::ID::PAR_2, ChemicalType::PROTEIN, Species::C_ELEGANS), 3.9e5);
    pInternalMedium->addMolecule(par2, float3(0, -1.f, 0));

    // Initialize maternal proteins at cell center
    float3 center(0.0f, 0.0f, 0.0f);

    // Add maternal CDK-1 and CYB-1 (Cyclin B)
    MPopulation cdk1(Molecule(StringDict::ID::CDK_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 1500.0);  // Initial amount above threshold (1000)
    MPopulation cyb1(Molecule(StringDict::ID::CYB_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 1500.0);  // Initial amount above threshold (1000)
    pInternalMedium->addMolecule(cdk1, center);
    pInternalMedium->addMolecule(cyb1, center);

    // Add centrosome-related proteins for proper centrosome function
    MPopulation cdk2(Molecule(StringDict::ID::CDK_2, ChemicalType::PROTEIN, Species::C_ELEGANS), 800.0);  // For centrosome duplication
    MPopulation cce1(Molecule(StringDict::ID::CCE_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 800.0);  // Cyclin E for centrosome duplication
    pInternalMedium->addMolecule(cdk2, center);
    pInternalMedium->addMolecule(cce1, center);
    
    // Add maternal ATP for translation
    pInternalMedium->addATP(50000.0, center);  // Sufficient ATP for early translation
    // Add maternal GTP/GDP nucleotide pools (optional explicit bookkeeping for GTPases)
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::GTP, ChemicalType::NUCLEOTIDE, Species::C_ELEGANS), 200000.0), center);
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::GDP, ChemicalType::NUCLEOTIDE, Species::C_ELEGANS), 200000.0), center);
    
    // Add maternal tRNAs (essential for translation bootstrap)
    // Without these, mRNAs (including tRNA mRNAs) cannot be translated
    addMaternalTRNAs(*pInternalMedium, center);

    // Maternal provisioning of polarity/contractility pathway components
    // Rho module: start mostly GDP-bound; allow dynamics to convert to GTP via ECT-2
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::RHO_1_GDP, ChemicalType::PROTEIN, Species::C_ELEGANS), 800000.0), center);
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::RHO_1_GTP, ChemicalType::PROTEIN, Species::C_ELEGANS), 200000.0), center);
    // ECT-2 (RhoGEF) and CHIN-1 (RhoGAP) as maternal proteins
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::ECT_2, ChemicalType::PROTEIN, Species::C_ELEGANS), 150000.0), center);
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::CHIN_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 120000.0), center);
    // CDC-42 module (initially GDP-biased)
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::CDC_42_GDP, ChemicalType::PROTEIN, Species::C_ELEGANS), 250000.0), center);
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::CDC_42_GTP, ChemicalType::PROTEIN, Species::C_ELEGANS), 50000.0), center);
    // Myosin II as contractility proxy (will be cortex-enriched by later mechanics)
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::NMY_2, ChemicalType::PROTEIN, Species::C_ELEGANS), 300000.0), center);
    // AIR-1 (Aurora A) maternally supplied; will enrich at centrosomes/MTs
    pInternalMedium->addMolecule(MPopulation(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 50000.0), center);

    return pInternalMedium;
}

Worm::Worm()
{
    // Initialize the string dictionary first
    StringDict::initialize();
    
    auto chromosomes = initializeGenes();
    auto pInternalMedium = createZygoteMedium();
    auto pCell = Cell::createCell(pInternalMedium, chromosomes, CellType::Zygote, Species::C_ELEGANS);
    
    // Simulate fertilization and seed maternal γ-tubulin near the posterior centrosome
    seedCentrosomeAndMaternalGammaTubulin(pCell, float3(0.0f, -0.8f, 0.0f));
    
    // Create CellSim to wrap the Cell
    auto pCellSim = std::make_shared<CellSim>(pCell);
    m_pCellSims.push_back(pCellSim);
    
    // Set up the data collector
    setupDataCollector();
}

void Worm::seedCentrosomeAndMaternalGammaTubulin(std::shared_ptr<Cell> pCell, const float3& posteriorEntryPoint)
{
    // Add centrosome at posterior entry point (sperm-derived centrioles)
    auto pCentrosome = std::make_shared<Centrosome>(std::weak_ptr<Cell>(pCell), posteriorEntryPoint);
    pCell->addOrganelle(StringDict::ID::ORGANELLE_CENTROSOME, pCentrosome);

    // Seed maternal γ-tubulin at the centrosome location to enable initial Y_TuRC formation pre-duplication
    auto& medium = pCell->getInternalMedium();
    MPopulation gammaTubulinSeed(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, Species::C_ELEGANS), 1000.0);
    medium.addMolecule(gammaTubulinSeed, posteriorEntryPoint);
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
    // Resolve data/simOutput via FileUtils and write sim.csv there
    std::filesystem::path simOutPath;
    std::string simCsv = "sim.csv";
    if (FileUtils::findTheFolder("data/simOutput", simOutPath)) {
        std::filesystem::create_directories(simOutPath);
        simCsv = (simOutPath / simCsv).string();
    }
    m_pDataCollector = std::make_unique<DataCollector>(
        internalMedium, 
        simCsv,
        5.0f  // Collect data every 5 seconds
    );
    // Provide cell for global metrics and enable nucleation site tracking
    m_pDataCollector->setCell(m_pCellSims[0]->getCell());
    m_pDataCollector->setTrackNucleationSites(true);
    
    // Add single collection point at posterior for γ-tubulin protein and mRNA
    float3 posteriorPos(0.0f, -0.8f, 0.0f); // Posterior region near centrosome
    m_pDataCollector->addCollectionPoint(
        posteriorPos,
        "Posterior",
        {
            Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, Species::C_ELEGANS),
            Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::MRNA, Species::C_ELEGANS)
        }
    );

    // Add cortex-bound PAR protein sampling for polarization analysis
    std::string par3Membrane = MoleculeWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::ORGANELLE_CORTEX);
    std::string par2Membrane = MoleculeWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::ORGANELLE_CORTEX);

    float3 anteriorCortex(0.0f, 0.9f, 0.0f);
    float3 posteriorCortex(0.0f, -0.9f, 0.0f);

    m_pDataCollector->addCollectionPoint(
        anteriorCortex,
        "AnteriorCortex",
        {
            Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, Species::C_ELEGANS),
            Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, Species::C_ELEGANS)
        }
    );

    m_pDataCollector->addCollectionPoint(
        posteriorCortex,
        "PosteriorCortex",
        {
            Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, Species::C_ELEGANS),
            Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, Species::C_ELEGANS)
        }
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
    
    // Interval-based data collection (every configured seconds)
    if (m_pDataCollector) {
        m_pDataCollector->update(m_fTotalTime);
    }
}

bool Worm::validatePARPolarization(float fTimeSec) const
{
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    
    // Sample at cortex-aligned positions
    float3 anteriorPos(0.0f, 0.9f, 0.0f);
    float3 posteriorPos(0.0f, -0.9f, 0.0f);
    
    // Get membrane-bound protein names using the utility function
    std::string par3Membrane = MoleculeWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::ORGANELLE_CORTEX);
    std::string par2Membrane = MoleculeWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::ORGANELLE_CORTEX);
    
    // Check membrane-bound proteins
    Species species = Species::C_ELEGANS;
    double anteriorPAR3 = internalMedium.getMoleculeNumber(Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, species), anteriorPos);
    double posteriorPAR3 = internalMedium.getMoleculeNumber(Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, species), posteriorPos);
    double anteriorPAR2 = internalMedium.getMoleculeNumber(Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, species), anteriorPos);
    double posteriorPAR2 = internalMedium.getMoleculeNumber(Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, species), posteriorPos);
    
    // Check during polarity establishment with ramped thresholds and robust ratios
    if (fTimeSec < POLARITY_ESTABLISHMENT_END_SEC) {
        const double eps = 1e-6;
        double par3Ratio = (anteriorPAR3 + eps) / (posteriorPAR3 + eps);
        double par2Ratio = (posteriorPAR2 + eps) / (anteriorPAR2 + eps);

        // Grace period and threshold ramp: no strict check before 60s, moderate by 120s, full by 360s
        if (fTimeSec < 60.0f) {
            return true;
        }

        double requiredRatio = 1.5;
        if (fTimeSec >= 180.0f && fTimeSec < POLARITY_ESTABLISHMENT_END_SEC) {
            // Linear ramp from 1.5 at 180s to 3.0 at 360s
            double t0 = 180.0;
            double t1 = static_cast<double>(POLARITY_ESTABLISHMENT_END_SEC);
            double alpha = (static_cast<double>(fTimeSec) - t0) / (t1 - t0);
            if (alpha < 0.0) alpha = 0.0;
            if (alpha > 1.0) alpha = 1.0;
            requiredRatio = 1.5 + alpha * (3.0 - 1.5);
        } else if (fTimeSec >= POLARITY_ESTABLISHMENT_END_SEC) {
            requiredRatio = 3.0;
        }

        if (par3Ratio < requiredRatio) {
            LOG_INFO("Warning: Insufficient anterior %s polarization (ratio %.2lf < %.2lf) at %.2lf sec", par3Membrane.c_str(), par3Ratio, requiredRatio, fTimeSec);
            return false;
        }

        if (par2Ratio < requiredRatio) {
            LOG_INFO("Warning: Insufficient posterior %s polarization (ratio %.2lf < %.2lf) at %.2lf sec", par2Membrane.c_str(), par2Ratio, requiredRatio, fTimeSec);
            return false;
        }
    }
    
    return true;
}

bool Worm::validateCellCycle(float fTimeSec) const
{
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    float3 nuclearPos(0.0f, 0.0f, 0.0f);
    double cdk1Level = internalMedium.getMoleculeNumber(Molecule(StringDict::ID::CDK_1, ChemicalType::PROTEIN), nuclearPos);

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

	// Biologically grounded guard: centriole duplication should not occur too early.
	// This check is intentionally conservative: we prevent duplication before 6 minutes.
	// Literature context (see data/prompts/mtLiterature.txt):
	// - Duplication is S-phase–restricted and occurs well after meiotic exit; an earliest realistic
	//   window at 20–22 °C is ≳10–12 min post-fertilization, with PCM/γ-tubulin maturation rising
	//   toward NEBD (~12–15 min) and metaphase ~15 min.
	//   (Sonneville et al. 2012, J Cell Biol; PMC3265957. Baumgart et al. 2019, J Cell Biol; PMID:31636117.)
	// - Cyclin E/CDK-2 is required for centrosome assembly and couples duplication competence to
	//   the cell cycle (Cowan & Hyman 2006, Nat Cell Biol; PMID:17115027).
	// - SPD-2 (CEP192) functions upstream to enable duplication/PCM maturation and γ-tubulin recruitment
	//   (Kemp et al. 2004, Dev Cell; PMID:15068791). Pathway: SPD-2 → ZYG-1 → SAS-6/5/4.
	// - Centrosome size/nucleation capacity scales with a limiting maternal PCM pool (Decker et al. 2011,
	//   Curr Biol; PMID:21802300), further arguing against fixed early-time duplication.
	// Note: We keep the 6-minute guard as a conservative lower bound until an explicit S-phase + markers
	// gate is implemented. Adjust with temperature if using time guards across conditions.
	static constexpr float MIN_DUPLICATION_TIME_SEC = 360.0f; // 6 minutes (conservative lower bound)
	if (fTimeSec < MIN_DUPLICATION_TIME_SEC && pCentrosome->isDuplicated()) {
		LOG_INFO("Warning: Centrosome duplicated too early at %.2lf sec (before %.2f sec)", fTimeSec, (double)MIN_DUPLICATION_TIME_SEC);
		return false;
	}
    
    // Check centrosome position during different phases
    float3 centrosomePos = pCentrosome->getNormalizedPosition();
    CellCycleState cellCycleState = m_pCellSims[0]->getCell()->getCellCycleState();
    
    switch (cellCycleState) {
        case CellCycleState::INTERPHASE:
            // Early interphase: allow posterior localization; enforce proximity later
            if (fTimeSec >= 180.0f) { // start enforcing after 3 minutes
                if (std::abs(centrosomePos.x) > 0.2f || std::abs(centrosomePos.y) > 0.2f || std::abs(centrosomePos.z) > 0.2f) {
                    LOG_INFO("Warning: Centrosome too far from nucleus during interphase at %.2lf sec", fTimeSec);
                    return false;
                }
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

bool Worm::validateGammaTubulinLevels(float fTimeSec) const
{
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    // Sample only at the expected centrosome position
    float3 posteriorCentrosome(0.0f, -0.8f, 0.0f);

    double gammaProtCentro = internalMedium.getMoleculeNumber(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN), posteriorCentrosome);
    double gammaMRNACentro = internalMedium.getMoleculeNumber(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::MRNA), posteriorCentrosome);

    // Simple early-time expectations: within first 5 min, protein remains very low; mRNA ramps up
    // Thresholds chosen conservatively for current model scaling
    if (fTimeSec < 60.0f) {
        if (gammaProtCentro > 5.0) {
            LOG_INFO("Warning: γ-tubulin protein unusually high early (%.2lf) at %.2lf sec", gammaProtCentro, fTimeSec);
            return false;
        }
        // mRNA should begin to appear by ~5-10s and exceed ~0.05 by 60s in current kinetics
        if (fTimeSec > 10.0f && gammaMRNACentro < 0.01) {
            LOG_INFO("Warning: γ-tubulin mRNA too low (%.5lf) at %.2lf sec", gammaMRNACentro, fTimeSec);
            return false;
        }
    }

    // Later expectations: by 6-10 min, centrosome protein should reach tens-hundreds range
    if (fTimeSec >= 360.0f) {
        if (gammaProtCentro < 50.0) {
            LOG_INFO("Warning: γ-tubulin protein low at centrosome (%.2lf) at %.2lf sec", gammaProtCentro, fTimeSec);
            return false;
        }
    }

    return true;
}
