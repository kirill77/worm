#pragma once

#include "Organelle.h"
#include "Medium.h"
#include "geometry/vectors/vector.h"
#include <random>

class Mitochondrion : public Organelle
{
private:
    double m_fNumber;  // Number of mitochondria in the population
    
    static constexpr double ATP_PRODUCTION_RATE = 1000;
    static constexpr double N_INITIAL_MITOCHONDRIA = 200;
    
    // Random number generator for ATP distribution
    std::mt19937 m_rng;

    // Generate ATP at n random positions in the medium
    void generateATP(Medium& medium, uint32_t n, double amount);
    // Helper method to generate a random position within the medium
    float3 generateRandomPosition();

public:
    Mitochondrion(std::weak_ptr<Cell> pCell);

    void update(double dt, Cell& cell) override;
};

