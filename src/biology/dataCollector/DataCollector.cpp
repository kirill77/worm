#include "pch.h"
#include "DataCollector.h"
#include "utils/csvFile/CSVFileWriter.h"
// Needed for nucleation site counting
#include "biology/organelles/Cell.h"
#include "biology/organelles/Centrosome.h"
#include "biology/organelles/Nucleus.h"
#include "biology/organelles/Y_TuRC.h"
#include "biology/organelles/Cortex.h"
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
    
    // Average MT length (µm)
    dataRow.push_back(computeAverageMTLength());
    // Approximate global α/β tubulin protein concentrations (cell-average)
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, Species::C_ELEGANS)));
    dataRow.push_back(computeApproxCellAverageConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, Species::C_ELEGANS)));

    // Microtubule tip diagnostics (averaged across active MTs)
    double avgTipTubDimer = 0.0;
    double avgVGrowEff = 0.0;
    double avgPCatEff = 0.0;
    double avgPResEff = 0.0;
    double avgDistToCortex = 0.0;
    double fracContact = 0.0;
    size_t mtCountDiag = 0;
    if (auto cellPtr = m_cell.lock()) {
        auto pCentro = std::dynamic_pointer_cast<Centrosome>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_CENTROSOME));
        auto pCortex = std::dynamic_pointer_cast<Cortex>(cellPtr->getOrganelle(StringDict::ID::ORGANELLE_CORTEX));
        if (pCentro && pCortex) {
            const auto& rings = pCentro->getRingComplexes();
            // Centrosome world position
            float3 centroWorld = pCortex->normalizedToWorld(pCentro->getNormalizedPosition());
            for (const auto& r : rings) {
                if (!r || !r->hasActiveMT()) continue;
                const float length = r->getMTLengthMicroM();
                const float3 tipWorld = centroWorld + r->getTipPosition();
                // Convert tip to normalized for sampling
                const float3 tipNorm = pCortex->worldToNormalized(tipWorld);
                // Sample local tubulin and AIR-1
                double tubAlpha = m_medium.getMoleculeConcentration(Molecule(StringDict::ID::ALPHA_TUBULIN, ChemicalType::PROTEIN, cellPtr->getSpecies()), tipNorm);
                double tubBeta  = m_medium.getMoleculeConcentration(Molecule(StringDict::ID::BETA_TUBULIN,  ChemicalType::PROTEIN, cellPtr->getSpecies()), tipNorm);
                double tubDimer = std::min(tubAlpha, tubBeta);
                double air1     = m_medium.getMoleculeConcentration(Molecule(StringDict::ID::AIR_1, ChemicalType::PROTEIN, cellPtr->getSpecies()), tipNorm);
                // Distance to cortex along the ray from centrosome to tip
                float3 dirToTip = tipWorld - centroWorld;
                float distToTip = sqrtf(dirToTip.x * dirToTip.x + dirToTip.y * dirToTip.y + dirToTip.z * dirToTip.z);
                if (distToTip > 0.0f) {
                    dirToTip = float3(dirToTip.x / distToTip, dirToTip.y / distToTip, dirToTip.z / distToTip);
                }
                Cortex::CortexRay ray(centroWorld, dirToTip);
                bool hit = pCortex->findClosestIntersection(ray);
                float maxLen = hit ? ray.getDistance() : 0.0f;
                bool contact = (maxLen > 0.0f) && (distToTip >= maxLen - 1e-4f);
                double distToCortex = (maxLen > 0.0f) ? std::max(0.0, static_cast<double>(maxLen - distToTip)) : 0.0;
                // Recompute effective vGrow and probabilities to log (mirror Y_TuRC constants)
                const double vGrowMax = 0.45;   // µm/s
                const double K_tub    = 50.0;   // molecules/µm^3 proxy
                const double K_air    = 10.0;
                double vGrow = vGrowMax * (tubDimer / (K_tub + tubDimer));
                double pCatFree = contact ? 0.6 : 0.12; // s^-1 baseline modified at contact
                double f_air_cat = 1.0 / (1.0 + air1 / K_air);
                double f_tub_cat = (K_tub / (K_tub + std::max(tubDimer, 0.0)));
                double f_tub_res = 1.0 + (std::max(tubDimer, 0.0) / (K_tub + std::max(tubDimer, 0.0)));
                double pCatEff = pCatFree * f_air_cat * f_tub_cat;
                double pResEff = 0.05 * f_tub_res; // baseline 0.05 s^-1

                avgTipTubDimer += tubDimer;
                avgVGrowEff    += vGrow;
                avgPCatEff     += pCatEff;
                avgPResEff     += pResEff;
                avgDistToCortex+= distToCortex;
                fracContact    += contact ? 1.0 : 0.0;
                ++mtCountDiag;
            }
        }
    }
    if (mtCountDiag > 0) {
        double inv = 1.0 / static_cast<double>(mtCountDiag);
        avgTipTubDimer *= inv;
        avgVGrowEff    *= inv;
        avgPCatEff     *= inv;
        avgPResEff     *= inv;
        avgDistToCortex*= inv;
        fracContact    *= inv;
    }
    dataRow.push_back(avgTipTubDimer);
    dataRow.push_back(avgVGrowEff);
    dataRow.push_back(avgPCatEff);
    dataRow.push_back(avgPResEff);
    dataRow.push_back(avgDistToCortex);
    dataRow.push_back(fracContact);
    
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
    
    // Average MT length column (µm) and global α/β-tubulin protein cell-average
    headers.push_back("AverageMTLength(µm)");
    headers.push_back("ALPHA_TUBULIN[/µm^3]");
    headers.push_back("BETA_TUBULIN[/µm^3]");
    
    // Diagnostic MT dynamics columns
    headers.push_back("AvgTipTubulinDimer[/µm^3]");
    headers.push_back("AvgMT_vGrowEff(µm/s)");
    headers.push_back("AvgMT_pCatEff(1/s)");
    headers.push_back("AvgMT_pResEff(1/s)");
    headers.push_back("AvgMT_DistToCortex(µm)");
    headers.push_back("FracMT_ContactingCortex");
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