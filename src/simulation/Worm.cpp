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

    // Create and add anterior proteins
    ProteinPopulation par3;
    par3.m_sName = "PAR-3";
    par3.m_fNumber = 3.9e5;
    pMedium->addProtein(par3, float3(0, 0.95f, 0));

    ProteinPopulation par6;
    par6.m_sName = "PAR-6";
    par6.m_fNumber = 3.9e5;
    pMedium->addProtein(par6, float3(0, 0.95f, 0));

    ProteinPopulation pkc3;
    pkc3.m_sName = "PKC-3";
    pkc3.m_fNumber = 3.9e5;
    pMedium->addProtein(pkc3, float3(0, 0.95f, 0));

    // Create and add posterior proteins
    ProteinPopulation par1;
    par1.m_sName = "PAR-1";
    par1.m_fNumber = 3.9e5;
    pMedium->addProtein(par1, float3(0, -0.95f, 0));

    ProteinPopulation par2;
    par2.m_sName = "PAR-2";
    par2.m_fNumber = 3.9e5;
    pMedium->addProtein(par2, float3(0, -0.95f, 0));

    return pMedium;
}

Worm::Worm()
{
    initializeGenes();
    std::shared_ptr<Medium> pMedium = createZygoteMedium();
    auto pCell = std::make_shared<Cell>(pMedium);
    m_pCells.push_back(pCell);
}
