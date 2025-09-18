#include "GeneWiki.h"
#include "chemistry/molecules/StringDict.h"
#include "chemistry/molecules/Molecule.h"
#include "utils/fileUtils/fileUtils.h"
#include "utils/log/ILog.h"
#include "HttpClient.h"
#include "chemistry/molecules/SpeciesUtils.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <assert.h>

namespace {
    static std::string makeRepeatedCodonSequence(const char* codon, int repeats)
    {
        std::string s;
        s.reserve(static_cast<size_t>(repeats) * 3u);
        for (int i = 0; i < repeats; ++i)
        {
            s.push_back(codon[0]);
            s.push_back(codon[1]);
            s.push_back(codon[2]);
        }
        return s;
    }

    static bool getBuiltinTrnaSequence(StringDict::ID id, std::string& outSequence)
    {
        switch (id)
        {
        case StringDict::ID::TRNA_MET_ATG: outSequence = makeRepeatedCodonSequence("ATG", 25); return true;
        case StringDict::ID::TRNA_GLY_GGA: outSequence = makeRepeatedCodonSequence("GGA", 25); return true;
        case StringDict::ID::TRNA_GLY_GGT: outSequence = makeRepeatedCodonSequence("GGT", 25); return true;
        case StringDict::ID::TRNA_ALA_GCA: outSequence = makeRepeatedCodonSequence("GCA", 25); return true;
        case StringDict::ID::TRNA_ALA_GCC: outSequence = makeRepeatedCodonSequence("GCC", 25); return true;
        case StringDict::ID::TRNA_LEU_CTG: outSequence = makeRepeatedCodonSequence("CTG", 25); return true;
        case StringDict::ID::TRNA_LEU_CTC: outSequence = makeRepeatedCodonSequence("CTC", 25); return true;
        case StringDict::ID::TRNA_SER_TCA: outSequence = makeRepeatedCodonSequence("TCA", 25); return true;
        case StringDict::ID::TRNA_SER_TCG: outSequence = makeRepeatedCodonSequence("TCG", 25); return true;
        case StringDict::ID::TRNA_VAL_GTG: outSequence = makeRepeatedCodonSequence("GTG", 25); return true;
        case StringDict::ID::TRNA_VAL_GTC: outSequence = makeRepeatedCodonSequence("GTC", 25); return true;
        case StringDict::ID::TRNA_PRO_CCA: outSequence = makeRepeatedCodonSequence("CCA", 25); return true;
        case StringDict::ID::TRNA_THR_ACA: outSequence = makeRepeatedCodonSequence("ACA", 25); return true;
        case StringDict::ID::TRNA_ASP_GAC: outSequence = makeRepeatedCodonSequence("GAC", 25); return true;
        case StringDict::ID::TRNA_GLU_GAG: outSequence = makeRepeatedCodonSequence("GAG", 25); return true;
        case StringDict::ID::TRNA_LYS_AAG: outSequence = makeRepeatedCodonSequence("AAG", 25); return true;
        case StringDict::ID::TRNA_ARG_CGA: outSequence = makeRepeatedCodonSequence("CGA", 25); return true;
        case StringDict::ID::TRNA_HIS_CAC: outSequence = makeRepeatedCodonSequence("CAC", 25); return true;
        case StringDict::ID::TRNA_PHE_TTC: outSequence = makeRepeatedCodonSequence("TTC", 25); return true;
        case StringDict::ID::TRNA_TYR_TAC: outSequence = makeRepeatedCodonSequence("TAC", 25); return true;
        case StringDict::ID::TRNA_CYS_TGC: outSequence = makeRepeatedCodonSequence("TGC", 25); return true;
        case StringDict::ID::TRNA_TRP_TGG: outSequence = makeRepeatedCodonSequence("TGG", 25); return true;
        case StringDict::ID::TRNA_ASN_AAC: outSequence = makeRepeatedCodonSequence("AAC", 25); return true;
        case StringDict::ID::TRNA_GLN_CAG: outSequence = makeRepeatedCodonSequence("CAG", 25); return true;
        case StringDict::ID::TRNA_ILE_ATC: outSequence = makeRepeatedCodonSequence("ATC", 25); return true;
        default: break;
        }
        return false;
    }
}

GeneWiki::GeneWiki()
{
    // Load caches for all known species
    for (int si = 0; si < static_cast<int>(Species::COUNT); ++si)
        loadMissingCache(static_cast<Species>(si));
    // Initialize built-in aliases for common non-canonical names
    // tRNAs: map generic names to species-specific lookup strings as needed
    // Protein genes: map vertebrate-style symbols to C. elegans symbols
    // Note: keys are mRNA molecules (species-aware)
    m_lookupAliases[Molecule(StringDict::ID::CCE_1, ChemicalType::MRNA, Species::C_ELEGANS)] = "cye-1";
    m_lookupAliases[Molecule(StringDict::ID::PLK_4, ChemicalType::MRNA, Species::C_ELEGANS)] = "zyg-1";
    m_lookupAliases[Molecule(StringDict::ID::GAMMA_TUBULIN, ChemicalType::MRNA, Species::C_ELEGANS)] = "tbg-1";
    m_lookupAliases[Molecule(StringDict::ID::NINEIN, ChemicalType::MRNA, Species::C_ELEGANS)] = "noca-1";
    // PERICENTRIN has no direct worm ortholog; prefer spd-2 or spd-5 depending on intent
    // Default to spd-2 as a scaffold component for sequence retrieval
    m_lookupAliases[Molecule(StringDict::ID::PERICENTRIN, ChemicalType::MRNA, Species::C_ELEGANS)] = "spd-2";

    // Human-specific aliases removed; use GENERIC for cross-species where applicable
}

GeneWiki& GeneWiki::getInstance()
{
    static GeneWiki instance;
    return instance;
}

const std::vector<std::pair<Molecule, uint32_t>>& GeneWiki::getGeneData(const Molecule& geneMolecule) const
{
    const Species species = geneMolecule.getSpecies();
    const std::string& geneName = geneMolecule.getName();
    assert(geneMolecule.getType() == ChemicalType::MRNA && "getGeneData expects an mRNA molecule");
    if (!ensureGeneDataComputed(geneMolecule))
    {
        throw std::runtime_error("GeneData could not be computed for: " + geneName);
    }
    return m_geneData.at(geneMolecule).m_trnaRequirements;
}

bool GeneWiki::hasGeneData(const Molecule& geneMolecule) const
{
    const Species species = geneMolecule.getSpecies();
    const std::string& geneName = geneMolecule.getName();
    assert(geneMolecule.getType() == ChemicalType::MRNA && "hasGeneData expects an mRNA molecule");
    if (m_geneData.find(geneMolecule) != m_geneData.end())
        return true;
    return ensureGeneDataComputed(geneMolecule);
}

// Removed default sequences; sequences are loaded lazily from disk or fetched

std::filesystem::path GeneWiki::getGenesFolder() const
{
    std::filesystem::path genesPath;
    if (!FileUtils::findTheFolder("data/genes", genesPath))
    {
        // fallback to current working directory
        genesPath = std::filesystem::current_path() / "data" / "genes";
    }
    if (!std::filesystem::exists(genesPath))
    {
        std::error_code ec;
        std::filesystem::create_directories(genesPath, ec);
        if (ec)
        {
            LOG_WARN("Failed to create genes folder: %s", genesPath.string().c_str());
        }
    }
    return genesPath;
}

std::filesystem::path GeneWiki::getSpeciesFolder(Species species) const
{
    const std::filesystem::path base = getGenesFolder();
    // Use canonical species string for folder naming
    std::filesystem::path sf = base / std::filesystem::path(speciesToWString(species));
    if (!std::filesystem::exists(sf))
    {
        std::error_code ec;
        std::filesystem::create_directories(sf, ec);
        if (ec)
        {
            LOG_WARN("Failed to create species folder: %s", sf.string().c_str());
        }
    }
    return sf;
}

std::string GeneWiki::sanitizeGeneNameForFile(const std::string& geneName)
{
    std::string out;
    out.reserve(geneName.size());
    for (char c : geneName)
    {
        if ((c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || (c >= '0' && c <= '9') || c == '-')
            out.push_back(c);
        else
            out.push_back('_');
    }
    return out;
}

std::filesystem::path GeneWiki::getGeneFilePath(Species species, const std::string& geneName) const
{
    const std::filesystem::path folder = getSpeciesFolder(species);
    const std::string fileBase = sanitizeGeneNameForFile(geneName);
    // Species is encoded in the parent folder; filename need not include species prefix
    return folder / (std::string("g_") + fileBase + ".fa");
}

std::filesystem::path GeneWiki::getMissingCacheFilePath(Species species) const
{
    // store alongside genes in species folder
    const std::filesystem::path folder = getSpeciesFolder(species);
    return folder / "missing_genes.cache";
}

bool GeneWiki::loadSequenceFromFile(const std::filesystem::path& filePath, std::string& outSequence) const
{
    if (!std::filesystem::exists(filePath))
        return false;
    std::ifstream in(filePath, std::ios::in);
    if (!in.is_open())
        return false;
    std::stringstream ss;
    std::string line;
    while (std::getline(in, line))
    {
        if (!line.empty() && line[0] == '>')
            continue; // skip FASTA header
        for (char c : line)
        {
            if (c != '\r' && c != '\n' && c != ' ' && c != '\t')
                ss << (char)toupper(c);
        }
    }
    outSequence = ss.str();
    return !outSequence.empty();
}

bool GeneWiki::saveSequenceToFile(const std::filesystem::path& filePath, const std::string& sequence) const
{
    std::ofstream out(filePath, std::ios::out | std::ios::trunc);
    if (!out.is_open())
        return false;
    out << "> autogenerated; DO NOT EDIT BY HAND\n";
    // write as wrapped FASTA for readability
    const size_t wrap = 80;
    for (size_t i = 0; i < sequence.size(); i += wrap)
    {
        out << sequence.substr(i, std::min(wrap, sequence.size() - i)) << "\n";
    }
    return true;
}

std::string GeneWiki::resolveLookupName(const Molecule& mrna) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    auto it = m_lookupAliases.find(mrna);
    if (it != m_lookupAliases.end())
        return it->second;
    return mrna.getName();
}

bool GeneWiki::fetchSequenceFromPublicDb(const Molecule& mrna, std::string& outSequence) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    const Species species = mrna.getSpecies();
    const std::string lookupName = resolveLookupName(mrna);
    // Try Ensembl REST; species-specific endpoint selection
    // Currently implements C. elegans; fallback for Human can be added similarly
    // 1) Lookup to get stable ID
    std::wstring base = L"https://rest.ensembl.org";
    std::wstring speciesPath = speciesToWString(species);
    std::wstring lookup = base + L"/xrefs/symbol/" + speciesPath + L"/" + std::wstring(lookupName.begin(), lookupName.end()) + L"?content-type=application/json";
    auto r1 = HttpClient::get(lookup, { {L"User-Agent", L"worm/1.0"} });
    if (r1.statusCode != 200 || r1.body.empty())
    {
        LOG_WARN("Ensembl xrefs lookup failed (%d) for '%s'", r1.statusCode, lookupName.c_str());
    }
    // Parse minimal: find first "id":"..."
    std::string id;
    size_t posId = r1.body.find("\"id\":\"");
    if (posId != std::string::npos)
    {
        posId += 6;
        size_t end = r1.body.find('"', posId);
        if (end != std::string::npos)
            id = r1.body.substr(posId, end - posId);
    }
    if (!id.empty())
    {
        // 2) Fetch sequence
        std::wstring wid(id.begin(), id.end());
        std::wstring seqUrl = base + L"/sequence/id/" + wid + L"?content-type=text/x-fasta";
        auto r2 = HttpClient::get(seqUrl, { {L"User-Agent", L"worm/1.0"} });
        if (r2.statusCode == 200 && !r2.body.empty())
        {
            // Parse FASTA body to raw sequence
            std::stringstream ss(r2.body);
            std::string line, seq;
            while (std::getline(ss, line))
            {
                if (!line.empty() && line[0] == '>') continue;
                for (char c : line) if (c!='\r'&&c!='\n'&&c!=' '&&c!='\t') seq.push_back((char)toupper(c));
            }
            if (!seq.empty())
            {
                outSequence = std::move(seq);
                return true;
            }
        }
        LOG_WARN("Ensembl sequence fetch failed (%d) for '%s' id '%s'", r2.statusCode, lookupName.c_str(), id.c_str());
    }

    // No synthetic fallback: report failure so callers can skip creating interactions
    LOG_WARN("Sequence not found for %s gene '%s' in public DB; skipping.",
        speciesToString(species), lookupName.c_str());
    return false;
}

bool GeneWiki::loadSequence(const Molecule& mrna, std::string& outSequence) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    const Species species = mrna.getSpecies();
    const std::string& geneName = mrna.getName();

    // If previously marked missing, skip fetch attempts
    if (isMarkedMissing(mrna))
    {
        return false;
    }

    const std::filesystem::path p = getGeneFilePath(species, geneName);
    std::string seq;
    // Built-in tRNA sequences: short-circuit loading for tRNA genes
    if (getBuiltinTrnaSequence(mrna.getID(), seq))
    {
        outSequence = std::move(seq);
        markFound(mrna);
        return true;
    }
    if (loadSequenceFromFile(p, seq))
    {
        outSequence = std::move(seq);
        // ensure not marked missing anymore
        markFound(mrna);
        return true;
    }
    // Not on disk, attempt fetch then save
    if (!fetchSequenceFromPublicDb(mrna, seq))
    {
        // remember negative result
        markMissing(mrna);
        return false;
    }
    // ensure folder exists and save
    std::filesystem::create_directories(p.parent_path());
    if (!saveSequenceToFile(p, seq))
    {
        LOG_WARN("Failed to save gene file: %s", p.string().c_str());
    }
    outSequence = std::move(seq);
    markFound(mrna);
    return true;
}

StringDict::ID GeneWiki::codonToChargedTrnaId(const std::string &codon)
{
    if (codon == "ATG") return StringDict::ID::TRNA_MET_ATG_CHARGED;
    if (codon == "GGA") return StringDict::ID::TRNA_GLY_GGA_CHARGED;
    if (codon == "GGT") return StringDict::ID::TRNA_GLY_GGT_CHARGED;
    if (codon == "GCA") return StringDict::ID::TRNA_ALA_GCA_CHARGED;
    if (codon == "GCC") return StringDict::ID::TRNA_ALA_GCC_CHARGED;
    if (codon == "CTG") return StringDict::ID::TRNA_LEU_CTG_CHARGED;
    if (codon == "CTC") return StringDict::ID::TRNA_LEU_CTC_CHARGED;
    if (codon == "TCA") return StringDict::ID::TRNA_SER_TCA_CHARGED;
    if (codon == "TCG") return StringDict::ID::TRNA_SER_TCG_CHARGED;
    if (codon == "GTG") return StringDict::ID::TRNA_VAL_GTG_CHARGED;
    if (codon == "GTC") return StringDict::ID::TRNA_VAL_GTC_CHARGED;
    if (codon == "CCA") return StringDict::ID::TRNA_PRO_CCA_CHARGED;
    if (codon == "ACA") return StringDict::ID::TRNA_THR_ACA_CHARGED;
    if (codon == "GAC") return StringDict::ID::TRNA_ASP_GAC_CHARGED;
    if (codon == "GAG") return StringDict::ID::TRNA_GLU_GAG_CHARGED;
    if (codon == "AAG") return StringDict::ID::TRNA_LYS_AAG_CHARGED;
    if (codon == "CGA") return StringDict::ID::TRNA_ARG_CGA_CHARGED;
    if (codon == "CAC") return StringDict::ID::TRNA_HIS_CAC_CHARGED;
    if (codon == "TTC") return StringDict::ID::TRNA_PHE_TTC_CHARGED;
    if (codon == "TAC") return StringDict::ID::TRNA_TYR_TAC_CHARGED;
    if (codon == "TGC") return StringDict::ID::TRNA_CYS_TGC_CHARGED;
    if (codon == "TGG") return StringDict::ID::TRNA_TRP_TGG_CHARGED;
    if (codon == "AAC") return StringDict::ID::TRNA_ASN_AAC_CHARGED;
    if (codon == "CAG") return StringDict::ID::TRNA_GLN_CAG_CHARGED;
    if (codon == "ATC") return StringDict::ID::TRNA_ILE_ATC_CHARGED;
    return StringDict::ID::eUNKNOWN;
}

bool GeneWiki::ensureGeneDataComputed(const Molecule& mrna) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    const Molecule keyMol = mrna;
    if (m_geneData.find(keyMol) != m_geneData.end())
        return true;
    std::string seq;
    if (!loadSequence(keyMol, seq))
        return false;
    std::map<StringDict::ID, uint32_t> trnaCounts;
    for (size_t i = 0; i + 2 < seq.size(); i += 3)
    {
        std::string codon = seq.substr(i, 3);
        for (char &c : codon) c = (char)toupper(c);
        StringDict::ID trnaId = codonToChargedTrnaId(codon);
        if (trnaId == StringDict::ID::eUNKNOWN)
            continue;
        trnaCounts[trnaId] += 1u;
    }
    GeneData data;
    data.m_sequence = seq;
    data.m_trnaRequirements.reserve(trnaCounts.size());
    for (const auto &kv : trnaCounts)
    {
        Molecule trna(kv.first, ChemicalType::TRNA);
        data.m_trnaRequirements.emplace_back(trna, kv.second);
    }
    m_geneData[keyMol] = std::move(data);
    return true;
}

// makeKey removed: sequences are keyed by Molecule

void GeneWiki::loadMissingCache(Species species) const
{
    // Load a per-species cache file; file contains plain gene names
    const std::filesystem::path cachePath = getMissingCacheFilePath(species);
    if (!std::filesystem::exists(cachePath))
        return;
    std::ifstream in(cachePath);
    if (!in.is_open())
        return;
    std::string line;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        // Convert gene name to Molecule key (MRNA for this species)
        StringDict::ID id = StringDict::stringToId(line);
        Molecule keyMol(id, ChemicalType::MRNA, species);
        m_missingSequenceKeys.insert(keyMol);
    }
}

void GeneWiki::saveMissingCache(Species species) const
{
    // Rewrite the species file from the unified set; write plain gene names (no species prefix)
    const std::filesystem::path cachePath = getMissingCacheFilePath(species);
    std::ofstream out(cachePath, std::ios::out | std::ios::trunc);
    if (!out.is_open())
        return;
    for (const auto &keyMol : m_missingSequenceKeys)
    {
        if (keyMol.getSpecies() != species)
            continue;
        if (keyMol.getType() != ChemicalType::MRNA)
            continue;
        const std::string &gene = StringDict::idToString(keyMol.getID());
        out << gene << '\n';
    }
}

bool GeneWiki::isMarkedMissing(const Molecule& mrna) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    return m_missingSequenceKeys.find(mrna) != m_missingSequenceKeys.end();
}

void GeneWiki::markMissing(const Molecule& mrna) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    if (m_missingSequenceKeys.insert(mrna).second)
    {
        saveMissingCache(mrna.getSpecies());
    }
}

void GeneWiki::markFound(const Molecule& mrna) const
{
    assert(mrna.getType() == ChemicalType::MRNA);
    auto it = m_missingSequenceKeys.find(mrna);
    if (it != m_missingSequenceKeys.end())
    {
        m_missingSequenceKeys.erase(it);
        saveMissingCache(mrna.getSpecies());
    }
}