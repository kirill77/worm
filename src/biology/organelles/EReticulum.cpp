#include "pch.h"
#include "EReticulum.h"
#include "Cell.h"
#include "Medium.h"
#include <random>

EReticulum::EReticulum(std::weak_ptr<Cell> pCell)
    : Organelle(pCell)
{
}

void EReticulum::update(double dt, Cell& cell)
{
}

void EReticulum::synthesizeProteins(Medium& medium)
{
    // Add proteins to random positions in the cell
    float3 position = generateRandomPosition();
    MPopulation proteins("ER-Protein", PROTEIN_SYNTHESIS_RATE);
    medium.addProtein(proteins, position);
}

void EReticulum::synthesizeLipids(Medium& medium)
{
    // Add lipids to random positions in the cell
    float3 position = generateRandomPosition();
    MPopulation lipids("ER-Lipid", LIPID_SYNTHESIS_RATE);
    medium.addProtein(lipids, position);
}

float3 EReticulum::generateRandomPosition() const
{
    // Generate random position within the normalized range [-1, 1]
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<float> dist(-1.f, 1.f);
    
    return float3(dist(gen), dist(gen), dist(gen));
} 