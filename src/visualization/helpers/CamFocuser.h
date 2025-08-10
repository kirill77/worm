#pragma once

#include <memory>
#include <vector>
#include <string>
#include <cstdint>
#include "chemistry/StringDict.h"
#include "geometry/vectors/box.h"

// Forward declarations
class Organism;
class GPUCamera;
class CameraTransition;

/**
 * Handles camera focusing functionality for camera transitions to focus on specific objects
 */
class CamFocuser
{
public:
    /**
     * Constructor that initializes the focus list with default organelles
     */
    CamFocuser();

    /**
     * Creates a camera transition to focus on the next organelle in the focus list
     * 
     * @param pOrganism The organism containing cells to focus on
     * @param pCurrentCam The current camera to transition from
     * @param transitionDurationSec Duration of the camera transition in seconds
     * @return Shared pointer to CameraTransition, or nullptr if unable to create transition
     */
    std::shared_ptr<CameraTransition> goToNextFocus(
        std::shared_ptr<Organism> pOrganism,
        std::shared_ptr<GPUCamera> pCurrentCam,
        float transitionDurationSec = 1.0f);
    
    /**
     * Gets the display name of the organelle that was last focused on
     * 
     * @return Display name of the last focused organelle
     */
    std::string getLastFocusedOrganelleName() const;

private:
    /**
     * Gets the bounding box of a specific organelle in the first cell
     * 
     * @param pOrganism The organism containing the cells
     * @param organelleId The ID of the organelle to get the bounding box for
     * @return Bounding box of the organelle, or empty box if not found
     */
    box3 getBox(std::shared_ptr<Organism> pOrganism, StringDict::ID organelleId);

    // Focus management
    std::vector<StringDict::ID> m_focusList;  // List of organelles to cycle through
    uint32_t m_curFocus;                      // Current focus index in the list
};
