#pragma once

#include <vector>
#include <memory>
#include <string>
#include "MoleculeInteraction.h"

// Repository for molecule interaction data (separate from MoleculeWiki)
class InteractionsWiki
{
public:
	// Initialize interactions by loading from data files
	static void Initialize();

	// Get all known interactions
	static const std::vector<std::shared_ptr<MoleculeInteraction>>& GetMoleculeInteractions();

private:
	// Stored interactions
	static std::vector<std::shared_ptr<MoleculeInteraction>> s_moleculeInteractions;
};


