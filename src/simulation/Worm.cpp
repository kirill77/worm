#include "pch.h"
#include "Worm.h"
#include "Cell.h"
#include "Protein.h"
#include "Medium.h"

Worm::Worm()
{
    std::shared_ptr<Medium> pMedium = std::make_shared<Medium>();

    // add anterior proteins
    std::shared_ptr<Protein> pProtein = std::make_shared<Protein>();
    pProtein->m_sName = "PAR-3";
    pProtein->m_fNumber = 3.9e5;
    pMedium->addProtein(pProtein, float3(0, 0.95f, 0));

    pProtein = std::make_shared<Protein>();
    pProtein->m_sName = "PAR-6";
    pProtein->m_fNumber = 3.9e5;
    pMedium->addProtein(pProtein, float3(0, 0.95f, 0));

    pProtein = std::make_shared<Protein>();
    pProtein->m_sName = "PKC-3";
    pProtein->m_fNumber = 3.9e5;
    pMedium->addProtein(pProtein, float3(0, 0.95f, 0));

    // add posterior proteins
    pProtein = std::make_shared<Protein>();
    pProtein->m_sName = "PAR-1";
    pProtein->m_fNumber = 3.9e5;
    pMedium->addProtein(pProtein, float3(0, -0.95f, 0));

    pProtein = std::make_shared<Protein>();
    pProtein->m_sName = "PAR-2";
    pProtein->m_fNumber = 3.9e5;
    pMedium->addProtein(pProtein, float3(0, -0.95f, 0));

    auto pCell = std::make_shared<Cell>(pMedium);
    m_pCells.push_back(pCell);
}
