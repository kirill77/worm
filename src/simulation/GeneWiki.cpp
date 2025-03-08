#include "pch.h"
#include "GeneWiki.h"
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
    m_sequences["pie-1"] = "ATGCCGAATTCGTCGAATCCG";  // Germline specification
    m_sequences["pal-1"] = "ATGAATTCGCCGAATCCGTCG";  // Posterior fate
    m_sequences["skn-1"] = "ATGCCGTCGAATTCGAATCCG";  // Endoderm specification
    m_sequences["mex-3"] = "ATGTCGCCGAATTCGAATCCG";  // Anterior fate
    
    // Cell division and timing genes
    m_sequences["cdk-1"] = "ATGCCGAATTCGTCGAATCCG";  // Cell cycle control
    m_sequences["cyb-1"] = "ATGAATTCGCCGTCGAATCCG";  // Cyclin B
    m_sequences["plk-1"] = "ATGCCGTCGAATTCGAATCCG";  // Polo-like kinase
} 