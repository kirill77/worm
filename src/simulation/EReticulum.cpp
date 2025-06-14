#include "pch.h"
#include "EReticulum.h"
#include "Cell.h"
#include <random>

EReticulum::EReticulum(std::weak_ptr<Cell> pCell)
    : Organelle(pCell)
{
}

void EReticulum::update(double dt, Cell& cell)
{
    auto& internalMedium = cell.getInternalMedium();

    // Calculate synthesis amounts based on time step
    double proteinAmount = PROTEIN_SYNTHESIS_RATE * dt;
    double lipidAmount = LIPID_SYNTHESIS_RATE * dt;

    // Calculate ATP costs
    double proteinATPCost = proteinAmount * ATP_COST_PER_PROTEIN;
    double lipidATPCost = lipidAmount * ATP_COST_PER_LIPID;
    double totalATPCost = proteinATPCost + lipidATPCost;

    // Check if we have enough ATP
    if (cell.consumeATP(totalATPCost))
    {
        // Synthesize proteins and lipids
        synthesizeProteins(internalMedium);
        synthesizeLipids(internalMedium);
    }
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