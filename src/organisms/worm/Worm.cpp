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
#include "chemistry/molecules/simConstants.h"
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
    pDNA5->addGene(StringDict::ID::ALPHA_TUBULIN, MoleculeConstants::ALPHA_TUBULIN_EXPRESSION_RATE, MoleculeConstants::ALPHA_TUBULIN_BASAL_LEVEL); // α-tubulin (tba-1): cytoskeletal dimer component (1000x expression)
    pDNA5->addGene(StringDict::ID::BETA_TUBULIN,  MoleculeConstants::BETA_TUBULIN_EXPRESSION_RATE,  MoleculeConstants::BETA_TUBULIN_BASAL_LEVEL);  // β-tubulin (tbb-2): cytoskeletal dimer component (1000x expression)
    pDNA5->addGene(StringDict::ID::GAMMA_TUBULIN, 0.1, 0.05);  // γ-tubulin (tbg-1): nucleation scaffold
    
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
    // Endogenous tRNA production/export is active; no maternal provisioning needed.
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

    // Seed minimal PCM scaffold components to drive molecule-based maturation (biological placeholders)
    MPopulation spd2(Molecule(StringDict::ID::SPD_2, ChemicalType::PROTEIN, Species::C_ELEGANS), 300.0);
    MPopulation spd5(Molecule(StringDict::ID::SPD_5, ChemicalType::PROTEIN, Species::C_ELEGANS), 300.0);
    MPopulation plk1(Molecule(StringDict::ID::PLK_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 150.0);
    MPopulation air1(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, Species::C_ELEGANS), 100.0);
    medium.addMolecule(spd2, posteriorEntryPoint);
    medium.addMolecule(spd5, posteriorEntryPoint);
    medium.addMolecule(plk1, posteriorEntryPoint);
    medium.addMolecule(air1, posteriorEntryPoint);
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
    // Create a timestamp-based output folder under data/simOutput and write sim.csv there
    std::filesystem::path simOutPath;
    std::string simCsv = "sim.csv";
    if (FileUtils::getOrCreateSubFolderUsingTimestamp("data/simOutput", simOutPath)) {
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
    double anteriorPAR3 = internalMedium.getMoleculeConcentration(Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, species), anteriorPos);
    double posteriorPAR3 = internalMedium.getMoleculeConcentration(Molecule(StringDict::stringToId(par3Membrane), ChemicalType::PROTEIN, species), posteriorPos);
    double anteriorPAR2 = internalMedium.getMoleculeConcentration(Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, species), anteriorPos);
    double posteriorPAR2 = internalMedium.getMoleculeConcentration(Molecule(StringDict::stringToId(par2Membrane), ChemicalType::PROTEIN, species), posteriorPos);
    
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
    // Use concentration-based, relative validation: compare nuclear CDK-1 to a coarse cell-average
    const double cdk1Nuclear = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::CDK_1, ChemicalType::PROTEIN), nuclearPos);
    // Sample a few off-center points to approximate a cell-wide mean (grid-agnostic heuristic)
    const float s = 0.5f;
    const float3 samplePts[6] = { float3( s, 0, 0), float3(-s, 0, 0), float3(0,  s, 0), float3(0, -s, 0), float3(0, 0,  s), float3(0, 0, -s) };
    double cdk1Sum = 0.0;
    for (const auto& p : samplePts) {
        cdk1Sum += internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::CDK_1, ChemicalType::PROTEIN), p);
    }
    const double cdk1Mean = cdk1Sum / 6.0;
    const double epsConc = 1e-12;
    const double cdk1Ratio = (cdk1Mean > epsConc) ? (cdk1Nuclear / cdk1Mean) : 0.0;

    // Before nuclear envelope breakdown (0-12.5 minutes): CDK-1 nuclear enrichment should be modest
    // Use ratio thresholds (unitless) instead of absolute counts
    static constexpr double PRE_NEBD_MAX_RATIO = 1.5;   // tuneable
    static constexpr double ENTRY_MIN_RATIO   = 2.0;   // tuneable
    if (fTimeSec < NUCLEAR_ENVELOPE_BREAKDOWN_SEC && cdk1Ratio > PRE_NEBD_MAX_RATIO) {
        LOG_INFO("Warning: CDK-1 nuclear enrichment high before NEBD (ratio %.2lf > %.2lf) at %.2lf sec", cdk1Ratio, PRE_NEBD_MAX_RATIO, fTimeSec);
        return false;
    }
    
    // During mitotic entry (12.5-15 minutes): CDK-1 should increase (higher nuclear enrichment)
    if (fTimeSec >= NUCLEAR_ENVELOPE_BREAKDOWN_SEC && fTimeSec < SPINDLE_ASSEMBLY_START_SEC && cdk1Ratio < ENTRY_MIN_RATIO) {
        LOG_INFO("Warning: CDK-1 nuclear enrichment low during mitotic entry (ratio %.2lf < %.2lf) at %.2lf sec", cdk1Ratio, ENTRY_MIN_RATIO, fTimeSec);
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
    // Validation data context (experimental references):
    // - C. elegans one-cell embryo has no calibrated absolute γ-tubulin (TBG-1) copy numbers per centrosome
    //   reported to date; therefore we validate using concentration trends rather than absolute counts.
    //   Maternal provisioning yields early centrosomal TBG-1 protein while tbg-1 mRNA is broadly cytoplasmic.
    //   (DeMeyer & Song 2017, microPublication Biology; DOI: 10.17912/W2CW8H; PMID: 32550353)
    // - At metaphase, C. elegans zygote centrosomes concentrate tubulin strongly and nucleate >10,000 MTs per
    //   centrosome; local α/β-tubulin reaches ~470 µM soluble + ~230 µM polymer (~660 µM total), implying
    //   robust γ-tubulin/γ-TuRC presence during maturation. (Baumgart et al. 2019, J Cell Biol; PMID: 31636117)
    // - Cross-species quantitative bounds used qualitatively: human mitotic centrosomes contain ~1,340 γ-tubulin
    //   copies/centrosome (Bauer et al. 2016, EMBO J; PMID: 27539480), and interphase levels are ~5–20% of mitotic
    //   (Haren 2023 review, J Cell Biol; PMID: 37695451). These inform expectations for rising γ-tubulin near PCM
    //   later in the cycle but are not enforced as hard counts in C. elegans.
    // Policy derived from these data:
    // - Early (<60 s): allow relatively high γ-tubulin protein concentration due to maternal seeding; only flag
    //   extreme outliers. Require mRNA concentration to begin rising after ~10 s.
    // - Later (≥360 s): expect nonzero γ-tubulin protein concentration near the posterior centrosome as PCM matures.
    auto& internalMedium = m_pCellSims[0]->getCell()->getInternalMedium();
    // Sample only at the expected centrosome position
    float3 posteriorCentrosome(0.0f, -0.8f, 0.0f);

    double gammaProtCentro = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::PROTEIN, Species::C_ELEGANS), posteriorCentrosome);
    double gammaMRNACentro = internalMedium.getMoleculeConcentration(Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::MRNA, Species::C_ELEGANS), posteriorCentrosome);

    // Simple early-time expectations: within first 5 min, protein concentration remains low; mRNA concentration ramps up
    // Use concentration (molecules per µm^3); remove strict early protein cap (maternal provisioning)
    if (fTimeSec < 60.0f) {
        if (gammaProtCentro > 1e6) {
            LOG_INFO("Warning: γ-tubulin protein concentration extremely high early (%.8lf /µm^3) at %.2lf sec", gammaProtCentro, fTimeSec);
            return false;
        }
        // mRNA should begin to appear by ~5-10s and exceed a minimal concentration by 60s
        if (fTimeSec > 10.0f && gammaMRNACentro < 1e-6) {
            LOG_INFO("Warning: γ-tubulin mRNA concentration too low (%.8lf /µm^3) at %.2lf sec", gammaMRNACentro, fTimeSec);
            return false;
        }
    }

    // Later expectations: by 6-10 min, centrosome protein concentration should exceed a minimal threshold
    if (fTimeSec >= 360.0f) {
        if (gammaProtCentro < 1e-6) {
            LOG_INFO("Warning: γ-tubulin protein concentration low at centrosome (%.8lf /µm^3) at %.2lf sec", gammaProtCentro, fTimeSec);
            return false;
        }
    }

    return true;
}
