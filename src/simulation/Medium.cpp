#include "pch.h"
#include "Medium.h"


SubMedium& Medium::findBox(const float3& p)
{
    uint32_t index = 0;
    for (uint32_t u = 0; u < 3; ++u)
    {
        assert(p[u] >= -1 && p[u] <= 1);
        float fNormalized = (p[u] + 1) / 2;
        uint32_t tmp = std::min(GRID_RES - 1, (uint32_t)(GRID_RES * fNormalized));
        index = (index * GRID_RES) + tmp;
    }
    return m_pBoxes[index];
}

void Medium::addProtein(std::shared_ptr<Protein> pProtein, const float3& where)
{
    findBox(where).m_pProteins.push_back(pProtein);
}