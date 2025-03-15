#include "pch.h"
#include "Worm.h"
#include "simulation/Cell.h"
#include "simulation/Protein.h"
#include "simulation/Medium.h"
#include "simulation/Membrane.h"
#include "simulation/Spindle.h"
#include "log/ILog.h"

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
    pDNA1->addGene("mex-3", 0.8, 0.1);  // Anterior fate
    pDNA1->addGene("plk-1", 1.2, 0.2);  // Polo-like kinase

    // Chromosome II
    pDNA2->addGene("skn-1", 0.8, 0.1);  // Endoderm specification
    pDNA2->addGene("cyb-1", 1.2, 0.2);  // Cyclin B

    // Chromosome III
    pDNA3->addGene("pal-1", 0.8, 0.1);  // Posterior fate
    pDNA3->addGene("cdk-1", 1.2, 0.2);  // Cell cycle control

    // Chromosome IV
    pDNA4->addGene("pie-1", 0.8, 0.1);  // Germline specification

    // Create chromosomes with their respective DNA
    chromosomes.emplace_back(pDNA1);
    chromosomes.emplace_back(pDNA2);
    chromosomes.emplace_back(pDNA3);
    chromosomes.emplace_back(pDNA4);
    chromosomes.emplace_back(pDNA5);  // Empty for now
    chromosomes.emplace_back(pDNA6);  // Empty for now

    return chromosomes;
}

std::shared_ptr<Membrane> Worm::createZygoteMembrane()
{
    // Create the internal medium
    std::shared_ptr<Medium> pInternalMedium = std::make_shared<Medium>();

    // Create and add anterior proteins at the anterior cortex
    ProteinPopulation par3("PAR-3", 3.9e5);
    pInternalMedium->addProtein(par3, float3(0, 0.95f, 0));

    ProteinPopulation par6("PAR-6", 3.9e5);
    pInternalMedium->addProtein(par6, float3(0, 0.95f, 0));

    ProteinPopulation pkc3("PKC-3", 3.9e5);
    pInternalMedium->addProtein(pkc3, float3(0, 0.95f, 0));

    // Create and add posterior proteins at the posterior cortex
    ProteinPopulation par1("PAR-1", 3.9e5);
    pInternalMedium->addProtein(par1, float3(0, -0.95f, 0));

    ProteinPopulation par2("PAR-2", 3.9e5);
    pInternalMedium->addProtein(par2, float3(0, -0.95f, 0));

    // Initialize maternal proteins at cell center
    float3 center(0.0f, 0.0f, 0.0f);

    // Add maternal CDK-1 and CYB-1 (Cyclin B)
    ProteinPopulation cdk1("CDK-1", 1500.0);  // Initial amount above threshold (1000)
    ProteinPopulation cyb1("CYB-1", 1500.0);  // Initial amount above threshold (1000)
    pInternalMedium->addProtein(cdk1, center);
    pInternalMedium->addProtein(cyb1, center);

    // Create a membrane with the internal medium
    return std::make_shared<Membrane>(pInternalMedium);
}

Worm::Worm()
{
    auto chromosomes = initializeGenes();
    std::shared_ptr<Membrane> pMembrane = createZygoteMembrane();
    auto pCell = std::make_shared<Cell>(pMembrane, chromosomes);
    m_pCells.push_back(pCell);
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
    
    float3 anteriorPos(0.0f, 0.8f, 0.0f);
    float3 posteriorPos(0.0f, -0.8f, 0.0f);
    
    double anteriorPAR3 = internalMedium.getProteinNumber("PAR-3", anteriorPos);
    double posteriorPAR3 = internalMedium.getProteinNumber("PAR-3", posteriorPos);
    double anteriorPAR2 = internalMedium.getProteinNumber("PAR-2", anteriorPos);
    double posteriorPAR2 = internalMedium.getProteinNumber("PAR-2", posteriorPos);
    
    // Check during polarity establishment (0-6 minutes)
    if (fTimeSec < POLARITY_ESTABLISHMENT_END_SEC) {
        if (anteriorPAR3 / (posteriorPAR3 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient anterior PAR-3 polarization at %.2lf sec", fTimeSec);
            return false;
        }
        
        if (posteriorPAR2 / (anteriorPAR2 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient posterior PAR-2 polarization at %.2lf sec", fTimeSec);
            return false;
        }
    }
    
    return true;
}

bool Worm::validateCellCycle(float fTimeSec) const
{
    auto& internalMedium = m_pCells[0]->getInternalMedium();
    float3 nuclearPos(0.0f, 0.0f, 0.0f);
    double cdk1Level = internalMedium.getProteinNumber("CDK-1", nuclearPos);

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
    
    auto pSpindle = m_pCells[0]->getSpindle();
    float3 spindlePos = pSpindle->getPosition();
    
    if (spindlePos.y > -0.1f) {
        LOG_INFO("Warning: Spindle not properly positioned toward posterior at %.2lf sec", fTimeSec);
        return false;
    }
    
    return true;
}
