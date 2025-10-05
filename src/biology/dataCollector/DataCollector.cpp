#include "pch.h"
#include "DataCollector.h"
#include "utils/csvFile/CSVFileWriter.h"
// Needed for nucleation site counting
#include "biology/organelles/Cell.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/Nucleus.h"
#include "biology/organelles/Y_TuRC.h"
#include "chemistry/molecules/StringDict.h"

DataCollector::DataCollector(Medium& medium, const std::string& outputFile, double collectionInterval)
    : m_medium(medium)
    , m_outputFile(outputFile)
    , m_collectionInterval(collectionInterval)
{
    // Initialize performance metrics
    m_performanceMetrics["StepTime"] = 0.0;
    
    // Defer CSV file creation until we have collection points and know the headers
}

void DataCollector::addCollectionPoint(const float3& position, const std::string& name, 
                                     const std::vector<Molecule>& molecules)
{
    // Add new collection point
    CollectionPoint point;
    point.position = position;
    point.name = name;
    point.molecules = molecules;
    
    m_collectionPoints.push_back(point);
    
    // Create or recreate the CSV file with updated headers
    if (!m_csvFile) {
        // First time creation
        m_csvFile = std::make_shared<CSVFileWriter>(m_outputFile, generateHeaders());
        m_csvFile->setPrecision(6);
    } else {
        // Recreate with new headers
        std::string filename = m_outputFile;
        m_csvFile = std::make_shared<CSVFileWriter>(filename, generateHeaders());
        m_csvFile->setPrecision(6);
        // Note: This will overwrite previous data
    }
}

bool DataCollector::update(double currentTime)
{
    // Check if it's time to collect data
    if (currentTime >= m_lastCollectionTime + m_collectionInterval) {
        // Compute wall-clock delta since last collection
        double stepMs = 0.0;
        auto now = std::chrono::high_resolution_clock::now();
        if (m_hasLastWallTime) {
            auto micros = std::chrono::duration_cast<std::chrono::microseconds>(now - m_lastWallTime).count();
            stepMs = static_cast<double>(micros) / 1000.0;
        }
        m_lastWallTime = now;
        m_hasLastWallTime = true;

        collectData(currentTime, stepMs);
        return true;
    }
    return false;
}

void DataCollector::forceCollection(double currentTime, double stepTimeMs)
{
    collectData(currentTime, stepTimeMs);
}

void DataCollector::collectData(double currentTime, double stepTimeMs)
{
    // Make sure we have a CSV file and collection points
    if (!m_csvFile || m_collectionPoints.empty()) {
        return;
    }
    
    // Store the step time
    m_performanceMetrics["StepTime"] = stepTimeMs;
    
    // Prepare row of data
    std::vector<double> dataRow;
    dataRow.push_back(currentTime);
    if (m_trackNucleationSites) {
        dataRow.push_back(computeRingCount());
    }
    
    dataRow.push_back(computeAverageMTLength());
    // Global average alpha/beta tubulin protein concentrations (approximate)
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, Species::C_ELEGANS)));
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, Species::C_ELEGANS)));
    // Global average alpha/beta tubulin mRNA concentrations (approximate)
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::MRNA, Species::C_ELEGANS)));
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::MRNA, Species::C_ELEGANS)));
    // Centrosome-proximal alpha/beta tubulin (posterior centrosome)
    float3 posteriorCentro(0.0f, -0.8f, 0.0f);
    Species cel = Species::C_ELEGANS;
    dataRow.push_back(m_medium.getMoleculeConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, cel), posteriorCentro));
    dataRow.push_back(m_medium.getMoleculeConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, cel), posteriorCentro));
    dataRow.push_back(m_medium.getMoleculeConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::MRNA, cel), posteriorCentro));
    dataRow.push_back(m_medium.getMoleculeConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::MRNA, cel), posteriorCentro));
    
    // Probe tRNA: nuclear counts and cytosolic center concentrations
    double trnaUn_nuc = 0.0, trnaCh_nuc = 0.0;
    if (auto cellPtr = m_cell.lock()) {
        auto pNucleus = std::dynamic_pointer_cast<Nucleus>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_NUCLEUS));
        if (pNucleus) {
            trnaUn_nuc = pNucleus->getNuclearMoleculeAmount(m_probeTRNAUncharged);
            trnaCh_nuc = pNucleus->getNuclearMoleculeAmount(m_probeTRNACharged);
        }
    }
    dataRow.push_back(trnaUn_nuc);
    dataRow.push_back(trnaCh_nuc);
    float3 center(0,0,0);
    dataRow.push_back(m_medium.getMoleculeConcentration(m_probeTRNAUncharged, center));
    dataRow.push_back(m_medium.getMoleculeConcentration(m_probeTRNACharged, center));
    
    // Append step time per row (ms)
    dataRow.push_back(stepTimeMs);
    
    // Add the row to the CSV file
    m_csvFile->addRow(dataRow);
    
    // Update tracking variables
    m_lastCollectionTime = currentTime;
    m_dataPointCount++;
    
    // Flush every 10 data points to ensure data is written
    if (m_dataPointCount % 10 == 0) {
        m_csvFile->flush();
    }
}

std::vector<std::string> DataCollector::generateHeaders() const
{
    std::vector<std::string> headers;
    
    // First columns: time and, optionally, nucleation site count
    headers.push_back("Time(s)");
    if (m_trackNucleationSites) {
        headers.push_back("NucleationSites(Y_TuRC)");
    }
    
    // Average MT length column (µm)
    headers.push_back("AverageMTLength(µm)");
    headers.push_back("ALPHA_TUBULIN[/µm^3]");
    headers.push_back("BETA_TUBULIN[/µm^3]");
    headers.push_back("ALPHA_TUBULIN_mRNA[/µm^3]");
    headers.push_back("BETA_TUBULIN_mRNA[/µm^3]");
    headers.push_back("ALPHA_TUBULIN@Centro[/µm^3]");
    headers.push_back("BETA_TUBULIN@Centro[/µm^3]");
    headers.push_back("ALPHA_TUBULIN_mRNA@Centro[/µm^3]");
    headers.push_back("BETA_TUBULIN_mRNA@Centro[/µm^3]");
    
    // Append per-row step time header
    headers.push_back("TRNA_GLU_GAG_nuclear(count)");
    headers.push_back("TRNA_GLU_GAG_CHARGED_nuclear(count)");
    headers.push_back("TRNA_GLU_GAG_cytosol_center[/µm^3]");
    headers.push_back("TRNA_GLU_GAG_CHARGED_cytosol_center[/µm^3]");
    headers.push_back("StepTime(ms)");
    
    return headers;
} 

double DataCollector::computeAverageMTLength() const
{
    size_t mtCount = 0;
    double mtSum = 0.0;
    if (auto cellPtr = m_cell.lock()) {
        auto pCentro = std::dynamic_pointer_cast<Centrosome>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
        if (pCentro) {
            const auto& rings = pCentro->getRingComplexes();
            for (const auto& r : rings) {
                if (r && r->hasActiveMT()) {
                    mtSum += static_cast<double>(r->getMTLengthMicroM());
                    ++mtCount;
                }
            }
        }
    }
    return (mtCount > 0) ? (mtSum / static_cast<double>(mtCount)) : 0.0;
}

double DataCollector::computeRingCount() const
{
    double ringCount = 0.0;
    if (auto cellPtr = m_cell.lock()) {
        auto pCentro = std::dynamic_pointer_cast<Centrosome>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
        if (pCentro) {
            ringCount = static_cast<double>(pCentro->getRingComplexes().size());
        }
    }
    return ringCount;
}

double DataCollector::computeApproxCellAverageConcentration(const Molecule& molecule) const
{
    static const float3 pts[7] = { float3(0,0,0), float3(0.5f,0,0), float3(-0.5f,0,0), float3(0,0.5f,0), float3(0,-0.5f,0), float3(0,0,0.5f), float3(0,0,-0.5f) };
    double sum = 0.0;
    for (const auto& p : pts) {
        sum += m_medium.getMoleculeConcentration(molecule, p);
    }
    return sum / 7.0;
}