#include "pch.h"
#include "Worm.h"
#include "chemistry/StringDict.h"
#include "biology/Cell.h"
#include "biology/Centrosome.h"
#include "biology/CellTypes.h"
#include "chemistry/Molecule.h"
#include "biology/Medium.h"
#include "biology/Cortex.h"
#include "biology/Spindle.h"
#include "chemistry/ProteinWiki.h"
#include "utils/log/ILog.h"
#include <chrono> // For high_resolution_clock
#include <cmath>  // For std::abs

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
    // Chromosome I
    pDNA1->addGene(StringDict::ID::MEX_3, 0.8, 0.1);  // Anterior fate
    pDNA1->addGene(StringDict::ID::PLK_1, 1.2, 0.2);  // Polo-like kinase

    // Chromosome II
    pDNA2->addGene(StringDict::ID::SKN_1, 0.8, 0.1);  // Endoderm specification
    pDNA2->addGene(StringDict::ID::CYB_1, 1.2, 0.2);  // Cyclin B

    // Chromosome III
    pDNA3->addGene(StringDict::ID::PAL_1, 0.8, 0.1);  // Posterior fate
    pDNA3->addGene(StringDict::ID::CDK_1, 1.2, 0.2);  // Cell cycle control

    // Chromosome IV
    pDNA4->addGene(StringDict::ID::PIE_1, 0.8, 0.1);  // Germline specification

    // Create chromosomes with their respective DNA
    chromosomes.emplace_back(pDNA1);
    chromosomes.emplace_back(pDNA2);
    chromosomes.emplace_back(pDNA3);
    chromosomes.emplace_back(pDNA4);
    chromosomes.emplace_back(pDNA5);  // Empty for now
    chromosomes.emplace_back(pDNA6);  // Empty for now

    return chromosomes;
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
    auto pCentrosome = std::make_shared<Centrosome>(std::weak_ptr<Cell>(pCell), float3(0, 0, 0));
    pCell->addOrganelle(StringDict::ID::ORGANELLE_CENTROSOME, pCentrosome);
    
    m_pCells.push_back(pCell);
    
    // Set up the data collector
    setupDataCollector();
}

void Worm::setupDataCollector()
{
    // Get access to the cell's internal medium
    if (m_pCells.empty()) {
        LOG_ERROR("Cannot set up data collector: no cells available");
        return;
    }
    
    auto& internalMedium = m_pCells[0]->getInternalMedium();
    
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
    std::string par2Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::BS_CORTEX);
    std::string par3Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::BS_CORTEX);
    std::string bindingSiteCortex = StringDict::idToString(StringDict::ID::BS_CORTEX);
    
    // Add collection points with specific proteins to track
    m_pDataCollector->addCollectionPoint(
        anteriorPos, 
        "Anterior", 
        {par2Membrane, par3Membrane, StringDict::idToString(StringDict::ID::PAR_2), StringDict::idToString(StringDict::ID::PAR_3), StringDict::idToString(StringDict::ID::PKC_3), bindingSiteCortex}
    );
    
    m_pDataCollector->addCollectionPoint(
        posteriorPos, 
        "Posterior", 
        {par2Membrane, par3Membrane, StringDict::idToString(StringDict::ID::PAR_1), StringDict::idToString(StringDict::ID::PAR_2), bindingSiteCortex}
    );
    
    // Add collection point for centrosome tracking
    m_pDataCollector->addCollectionPoint(
        centerPos, 
        "Centrosome", 
        {StringDict::idToString(StringDict::ID::GAMMA_TUBULIN), StringDict::idToString(StringDict::ID::PERICENTRIN), StringDict::idToString(StringDict::ID::NINEIN), StringDict::idToString(StringDict::ID::PLK_4), StringDict::idToString(StringDict::ID::CDK_2), StringDict::idToString(StringDict::ID::CCE_1)}
    );
}

void Worm::simulateStep(double dt)
{
    // Use local variables for timing
    auto stepStartTime = std::chrono::high_resolution_clock::now();
    
    // Call the base class simulateStep to handle standard simulation
    Organism::simulateStep(dt);
    
    // Calculate time taken for this step using local variables
    auto stepEndTime = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stepEndTime - stepStartTime);
    double stepTimeMsec = duration.count() / 1000.0; // Convert to milliseconds
    
    // Update total simulation time
    m_fTotalTime += static_cast<float>(dt);
    
    // Force data collection at the current time (regardless of interval)
    if (m_pDataCollector) {
        m_pDataCollector->forceCollection(m_fTotalTime, stepTimeMsec);
    }
}

// Validation thresholds based on experimental data
static constexpr double ANTERIOR_POSTERIOR_RATIO_THRESHOLD = 3.0;
static constexpr double NUCLEAR_SIZE_THRESHOLD = 0.8;
static constexpr double ASYMMETRIC_DIVISION_RATIO = 0.6;

// Timing constants (in simulation steps, where each step is 0.1 seconds)
static constexpr uint32_t POLARITY_ESTABLISHMENT_END = 3600;    // 6 minutes
static constexpr uint32_t POLARITY_MAINTENANCE_END = 6000;      // 10 minutes
static constexpr uint32_t NUCLEAR_ENVELOPE_BREAKDOWN = 7500;    // 12.5 minutes
static constexpr uint32_t SPINDLE_ASSEMBLY_START = 9000;        // 15 minutes
static constexpr uint32_t DIVISION_START = 11000;               // 18.3 minutes

bool Worm::validatePARPolarization(float fTimeSec) const
{
    auto& internalMedium = m_pCells[0]->getInternalMedium();
    
    float3 anteriorPos(0.0f, 1.f, 0.0f);
    float3 posteriorPos(0.0f, -1.f, 0.0f);
    
    // Get membrane-bound protein names using the utility function
    std::string par3Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_3), StringDict::ID::BS_CORTEX);
    std::string par2Membrane = ProteinWiki::GetBoundProteinName(StringDict::idToString(StringDict::ID::PAR_2), StringDict::ID::BS_CORTEX);
    
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
    auto& internalMedium = m_pCells[0]->getInternalMedium();
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
    
    auto pSpindle = std::dynamic_pointer_cast<Spindle>(m_pCells[0]->getOrganelle(StringDict::ID::ORGANELLE_SPINDLE));
    float3 spindlePos = pSpindle->getPosition();
    
    if (spindlePos.y > -0.1f) {
        LOG_INFO("Warning: Spindle not properly positioned toward posterior at %.2lf sec", fTimeSec);
        return false;
    }
    
    return true;
}

bool Worm::validateCentrosomeBehavior(float fTimeSec) const
{
    auto pCentrosome = std::dynamic_pointer_cast<Centrosome>(m_pCells[0]->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
    if (!pCentrosome) {
        // Before fertilization, there should be no centrosome
        if (fTimeSec < 1.0f) {  // Assume fertilization happens within 1 second
            return true;  // This is expected before fertilization
        }
        LOG_INFO("Warning: No centrosome found in cell at %.2lf sec (after expected fertilization time)", fTimeSec);
        return false;
    }
    
    // Check centrosome position during different phases
    float3 centrosomePos = pCentrosome->getPosition();
    CellCycleState cellCycleState = m_pCells[0]->getCellCycleState();
    
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
