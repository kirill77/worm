#include "pch.h"
#include "Worm.h"
#include "Cell.h"
#include "Protein.h"
#include "Medium.h"
#include "DNA.h"

void Worm::initializeGenes()
{
    m_pDNA = std::make_shared<DNA>();

    // Cell fate specification genes
    m_pDNA->addGene("pie-1", 0.8, 0.1);  // Germline specification
    m_pDNA->addGene("pal-1", 0.8, 0.1);  // Posterior fate
    m_pDNA->addGene("skn-1", 0.8, 0.1);  // Endoderm specification
    m_pDNA->addGene("mex-3", 0.8, 0.1);  // Anterior fate
    
    // Cell division and timing genes
    m_pDNA->addGene("cdk-1", 1.2, 0.2);  // Cell cycle control
    m_pDNA->addGene("cyb-1", 1.2, 0.2);  // Cyclin B
    m_pDNA->addGene("plk-1", 1.2, 0.2);  // Polo-like kinase
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
    initializeGenes();
    std::shared_ptr<Medium> pMedium = createZygoteMedium();
    auto pCell = std::make_shared<Cell>(pMedium);
    m_pCells.push_back(pCell);
}
