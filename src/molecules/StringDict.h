#pragma once

#include <unordered_map>
#include <string>
#include <assert.h>
#include <vector>

struct StringDict
{
    enum class ID {
        eUNKNOWN,
        // Binding surface types
        BS_MEMBRANE,
        BS_CORTEX,
        BS_CENTROSOME,
        
        // PAR proteins (polarity establishment)
        PAR_1,
        PAR_2,
        PAR_3,
        PAR_6,
        PKC_3,
        
        // Cell cycle proteins
        CDK_1,
        CDK_2,
        CYB_1,
        CCE_1,
        PLK_1,
        PLK_4,
        
        // Centrosome proteins
        GAMMA_TUBULIN,
        PERICENTRIN,
        NINEIN,
        
        // Cell fate specification genes
        MEX_3,
        SKN_1,
        PAL_1,
        PIE_1,
        
        // TODO: add all IDs here
    };

    static void initialize();
    static const std::string& idToString(ID id)
    {
        assert(m_idToString.size() > 0);
        return m_idToString[static_cast<size_t>(id)];
    }
    static ID stringToId(const std::string& s)
    {
        assert(m_stringToId.size() > 0);
        auto it = m_stringToId.find(s);
        return (it == m_stringToId.end()) ? ID::eUNKNOWN : it->second;
    }

private:
    static std::vector<std::string> m_idToString;
    static std::unordered_map<std::string, ID> m_stringToId;
};

