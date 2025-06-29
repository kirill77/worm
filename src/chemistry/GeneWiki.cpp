#include "pch.h"
#include "GeneWiki.h"
#include "StringDict.h"
#include <stdexcept>

GeneWiki::GeneWiki()
{
    initializeDefaultSequences();
}

GeneWiki& GeneWiki::getInstance()
{
    static GeneWiki instance;
    return instance;
}

const std::string& GeneWiki::getSequence(const std::string& geneName) const
{
    auto it = m_sequences.find(geneName);
    if (it == m_sequences.end())
    {
        throw std::runtime_error("Gene sequence not found: " + geneName);
    }
    return it->second;
}

void GeneWiki::initializeDefaultSequences()
{
    // Cell fate specification genes
    m_sequences[StringDict::idToString(StringDict::ID::PIE_1)] = "ATGCCGAATTCGTCGAATCCG";  // Germline specification
    m_sequences[StringDict::idToString(StringDict::ID::PAL_1)] = "ATGAATTCGCCGAATCCGTCG";  // Posterior fate
    m_sequences[StringDict::idToString(StringDict::ID::SKN_1)] = "ATGCCGTCGAATTCGAATCCG";  // Endoderm specification
    m_sequences[StringDict::idToString(StringDict::ID::MEX_3)] = "ATGTCGCCGAATTCGAATCCG";  // Anterior fate
    
    // Cell division and timing genes
    m_sequences[StringDict::idToString(StringDict::ID::CDK_1)] = "ATGCCGAATTCGTCGAATCCG";  // Cell cycle control
    m_sequences[StringDict::idToString(StringDict::ID::CYB_1)] = "ATGAATTCGCCGTCGAATCCG";  // Cyclin B
    m_sequences[StringDict::idToString(StringDict::ID::PLK_1)] = "ATGCCGTCGAATTCGAATCCG";  // Polo-like kinase
} 