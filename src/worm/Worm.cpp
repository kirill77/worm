#include "pch.h"
#include "Worm.h"
#include "simulation/Cell.h"
#include "simulation/Protein.h"
#include "simulation/Medium.h"

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
