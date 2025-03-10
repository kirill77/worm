#include "pch.h"
#include "Worm.h"
#include "simulation/Cell.h"
#include "simulation/Protein.h"
#include "simulation/Medium.h"
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

std::shared_ptr<Medium> Worm::createZygoteMedium()
{
    std::shared_ptr<Medium> pMedium = std::make_shared<Medium>();

    // Create and add anterior proteins at the anterior cortex
    ProteinPopulation par3("PAR-3", 3.9e5);
    pMedium->addProtein(par3, float3(0, 0.95f, 0));

    ProteinPopulation par6("PAR-6", 3.9e5);
    pMedium->addProtein(par6, float3(0, 0.95f, 0));

    ProteinPopulation pkc3("PKC-3", 3.9e5);
    pMedium->addProtein(pkc3, float3(0, 0.95f, 0));

    // Create and add posterior proteins at the posterior cortex
    ProteinPopulation par1("PAR-1", 3.9e5);
    pMedium->addProtein(par1, float3(0, -0.95f, 0));

    ProteinPopulation par2("PAR-2", 3.9e5);
    pMedium->addProtein(par2, float3(0, -0.95f, 0));

    // Initialize maternal proteins at cell center
    float3 center(0.0f, 0.0f, 0.0f);

    // Add maternal CDK-1 and CYB-1 (Cyclin B)
    ProteinPopulation cdk1("CDK-1", 1500.0);  // Initial amount above threshold (1000)
    ProteinPopulation cyb1("CYB-1", 1500.0);  // Initial amount above threshold (1000)
    pMedium->addProtein(cdk1, center);
    pMedium->addProtein(cyb1, center);

    return pMedium;
}

Worm::Worm()
{
    auto chromosomes = initializeGenes();
    std::shared_ptr<Medium> pMedium = createZygoteMedium();
    auto pCell = std::make_shared<Cell>(pMedium, chromosomes);
    m_pCells.push_back(pCell);
}

bool Worm::validatePARPolarization(uint32_t timestep) const
{
    auto medium = m_pCells[0]->getMedium();
    
    // Check anterior PAR concentration
    float3 anteriorPos(0.0f, 0.8f, 0.0f);  // Near anterior cortex
    float3 posteriorPos(0.0f, -0.8f, 0.0f); // Near posterior cortex
    
    double anteriorPAR3 = medium->getProteinNumber("PAR-3", anteriorPos);
    double posteriorPAR3 = medium->getProteinNumber("PAR-3", posteriorPos);
    double anteriorPAR2 = medium->getProteinNumber("PAR-2", anteriorPos);
    double posteriorPAR2 = medium->getProteinNumber("PAR-2", posteriorPos);
    
    // Early polarization (0-200 timesteps)
    if (timestep < 200) {
        LOG_INFO("Timestep {} - PAR Polarization Check:", timestep);
        LOG_INFO("Anterior PAR-3: {}, Posterior PAR-3: {}", anteriorPAR3, posteriorPAR3);
        LOG_INFO("Anterior PAR-2: {}, Posterior PAR-2: {}", anteriorPAR2, posteriorPAR2);
        
        // PAR-3 should be higher in anterior
        if (anteriorPAR3 / (posteriorPAR3 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient anterior PAR-3 polarization");
            return false;
        }
        
        // PAR-2 should be higher in posterior
        if (posteriorPAR2 / (anteriorPAR2 + 1.0) < ANTERIOR_POSTERIOR_RATIO_THRESHOLD) {
            LOG_INFO("Warning: Insufficient posterior PAR-2 polarization");
            return false;
        }
    }
    
    return true;
}

bool Worm::validateCellCycle(uint32_t timestep) const
{
    auto medium = m_pCells[0]->getMedium();
    
    // Check CDK-1 levels
    float3 nuclearPos(0.0f, 0.0f, 0.0f);  // Assuming nucleus is roughly centered
    double cdk1Level = medium->getProteinNumber("CDK-1", nuclearPos);
    
    LOG_INFO("Timestep {} - Cell Cycle Check:", timestep);
    LOG_INFO("CDK-1 level: {}", cdk1Level);
    
    // Early interphase (0-300 timesteps): CDK-1 should be relatively low
    if (timestep < 300 && cdk1Level > 1000) {
        LOG_INFO("Warning: CDK-1 levels too high for early interphase");
        return false;
    }
    
    // M-phase entry (300-600 timesteps): CDK-1 should increase
    if (timestep >= 300 && timestep < 600 && cdk1Level < 1500) {
        LOG_INFO("Warning: CDK-1 levels too low for M-phase entry");
        return false;
    }
    
    return true;
}

bool Worm::validateAsymmetricDivision(uint32_t timestep) const
{
    // Only check during late stages
    if (timestep < 800) return true;
    
    auto pSpindle = m_pCells[0]->getSpindle();
    
    // Check spindle position (should be posterior-shifted)
    float3 spindlePos = pSpindle->getPosition();
    if (spindlePos.y > -0.1f) {  // Spindle should be posterior-shifted
        LOG_INFO("Warning: Spindle not properly positioned toward posterior");
        return false;
    }
    
    return true;
}
