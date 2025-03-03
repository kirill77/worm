#pragma once

#include <memory>
#include <array>
#include <vector>
#include "math/vector.h"

struct Protein;

struct SubMedium
{
    std::vector<std::shared_ptr<Protein>> m_pProteins;
};

// describes the medium where everything is going on
class Medium
{
    static const uint32_t GRID_RES = 3;
    SubMedium& findBox(const float3& p);
    std::array<SubMedium, GRID_RES * GRID_RES * GRID_RES> m_pBoxes;

public:
    void addProtein(std::shared_ptr<Protein> pProtein, const float3& where);
};

