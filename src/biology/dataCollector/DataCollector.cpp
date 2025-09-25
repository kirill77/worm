#include "pch.h"
#include "DataCollector.h"
#include "utils/csvFile/CSVFileWriter.h"
// Needed for nucleation site counting
#include "biology/organelles/Cell.h"
#include "biology/organelles/Centrosome.h"
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
    dataRow.push_back(currentTime); // First column is time
    if (m_trackNucleationSites) {
        // Count Y_TuRC instances from cell's centrosome(s)
        double ringCount = 0.0;
        if (auto cellPtr = m_cell.lock()) {
            auto pCentro = std::dynamic_pointer_cast<Centrosome>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
            if (pCentro) {
                ringCount = static_cast<double>(pCentro->getRingComplexes().size());
            }
        }
        dataRow.push_back(ringCount);
    }
    
    // For each collection point, get molecule concentrations
    for (const auto& point : m_collectionPoints) {
        // For each molecule at this point
        for (const auto& molecule : point.molecules) {
            double concentration = m_medium.getMoleculeNumber(molecule, point.position);
            dataRow.push_back(concentration);
        }
    }
    
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
    
    // For each collection point, add molecule headers
    for (const auto& point : m_collectionPoints) {
        for (const auto& molecule : point.molecules) {
            std::string moleculeTypeStr;
            switch (molecule.getType()) {
                case ChemicalType::PROTEIN: moleculeTypeStr = "PROT"; break;
                case ChemicalType::MRNA: moleculeTypeStr = "mRNA"; break;
                case ChemicalType::TRNA: moleculeTypeStr = "tRNA"; break;
                case ChemicalType::DNA: moleculeTypeStr = "DNA"; break;
                case ChemicalType::NUCLEOTIDE: moleculeTypeStr = "NUC"; break;
                default: moleculeTypeStr = "OTHER"; break;
            }
            headers.push_back(molecule.getName() + "(" + moleculeTypeStr + ")_" + point.name);
        }
    }
    
    // Append per-row step time header
    headers.push_back("StepTime(ms)");
    
    return headers;
} 