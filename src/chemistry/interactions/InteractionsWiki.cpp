#include "InteractionsWiki.h"
#include "chemistry/molecules/MoleculeWiki.h"
#include "MoleculeInteractionLoader.h"
#include "utils/log/ILog.h"
#include "utils/fileUtils/fileUtils.h"
#include <filesystem>

// Static storage
std::vector<std::shared_ptr<MoleculeInteraction>> InteractionsWiki::s_moleculeInteractions;

void InteractionsWiki::Initialize()
{
	// Clear existing interactions
	s_moleculeInteractions.clear();

	// Ensure molecule info is initialized before interactions are loaded
	MoleculeWiki::Initialize();

	// Discover data directory (mirror logic from MoleculeWiki)
	std::filesystem::path dataPath;
	bool dataFolderFound = false;

	if (std::filesystem::exists("data/proteinRules")) {
		dataPath = "data/proteinRules";
		dataFolderFound = true;
	}
	else if (FileUtils::findTheFolder("data", dataPath)) {
		dataPath /= "proteinRules";
		if (std::filesystem::exists(dataPath)) {
			dataFolderFound = true;
		}
	}

	if (!dataFolderFound) {
		std::vector<std::string> commonPaths = {
			"../data/proteinRules",
			"../../data/proteinRules",
			"../../../data/proteinRules"
		};
		for (const auto& path : commonPaths) {
			if (std::filesystem::exists(path)) {
				dataPath = path;
				dataFolderFound = true;
				break;
			}
		}
	}

	if (dataFolderFound) {
		LOG_INFO("Loading molecule interactions from %s", dataPath.string().c_str());
		s_moleculeInteractions = MoleculeInteractionLoader::LoadAllInteractions(dataPath.string());
		if (s_moleculeInteractions.empty()) {
			LOG_ERROR("No molecule interactions were loaded from CSV files.");
		}
	}
	else {
		LOG_ERROR("Interaction data directory not found. Using default hardcoded interactions.");
	}
}

const std::vector<std::shared_ptr<MoleculeInteraction>>& InteractionsWiki::GetMoleculeInteractions()
{
	return s_moleculeInteractions;
}

std::vector<std::shared_ptr<MoleculeInteraction>> InteractionsWiki::GetInteractionsByMechanism(MoleculeInteraction::Mechanism mechanism)
{
	std::vector<std::shared_ptr<MoleculeInteraction>> result;
	result.reserve(s_moleculeInteractions.size());
	for (const auto& interaction : s_moleculeInteractions) {
		if (interaction->getMechanism() == mechanism) {
			result.push_back(interaction);
		}
	}
	return result;
}


