#include "CamFocuser.h"
#include "biology/simulation/Organism.h"
#include "biology/simulation/CellSim.h"
#include "biology/organelles/Cell.h"
#include "biology/organelles/Organelle.h"
#include "visualization/gpu/GPUCamera.h"
#include "visualization/helpers/CameraTransition.h"
#include "visualization/gpu/IVisObject.h"
#include "chemistry/StringDict.h"
#include "geometry/vectors/box.h"

CamFocuser::CamFocuser()
    : m_curFocus(0)
{
    // Initialize the focus list with default organelles
    m_focusList.push_back(StringDict::ID::ORGANELLE_CENTROSOME);
    m_focusList.push_back(StringDict::ID::ORGANELLE_CORTEX);
}

std::shared_ptr<CameraTransition> CamFocuser::goToNextFocus(
    std::shared_ptr<Organism> pOrganism,
    std::shared_ptr<GPUCamera> pCurrentCam,
    float transitionDurationSec)
{
    // Check if focus list is empty
    if (m_focusList.empty())
        return nullptr;

    // Get the current organelle to focus on
    StringDict::ID currentOrganelleId = m_focusList[m_curFocus];
    
    // Advance to next focus (with wraparound)
    m_curFocus = (m_curFocus + 1) % m_focusList.size();

    // Get bounding box for current organelle
    box3 organelleBounds = getBox(pOrganism, currentOrganelleId);
    if (organelleBounds.isempty())
        return nullptr;

    // Validate current camera
    if (!pCurrentCam)
        return nullptr;

    // Create target camera using fitBoxToView
    auto targetCamera = std::make_shared<GPUCamera>(*pCurrentCam);
    if (!targetCamera->fitBoxToView(organelleBounds))
        return nullptr;

    // Create and return the camera transition
    return std::make_shared<CameraTransition>(pCurrentCam, targetCamera,
        transitionDurationSec, organelleBounds);
}

box3 CamFocuser::getBox(std::shared_ptr<Organism> pOrganism, StringDict::ID organelleId)
{
    // Guard: need organism and at least one cell
    if (!pOrganism)
        return box3::empty();

    const auto& cells = pOrganism->getCellSims();
    if (cells.empty())
        return box3::empty();

    auto pCell = cells[0]->getCell();
    if (!pCell)
        return box3::empty();

    // Get specified organelle
    auto pOrganelle = pCell->getOrganelle(organelleId);
    if (!pOrganelle)
        return box3::empty();

    // Get visualization object and mesh node bounding box
    auto pVisObject = pOrganelle->getVisObject();
    if (!pVisObject)
        return box3::empty();

    auto meshNode = pVisObject->getMeshNode();
    box3 organelleBounds = meshNode.getBoundingBox();
    
    return organelleBounds;
}

std::string CamFocuser::getLastFocusedOrganelleName() const
{
    if (m_focusList.empty())
        return "Unknown";
    
    // Get the previous focus index (the one that was just focused on)
    uint32_t lastFocusIndex = (m_curFocus == 0) ? static_cast<uint32_t>(m_focusList.size()) - 1 : m_curFocus - 1;
    StringDict::ID organelleId = m_focusList[lastFocusIndex];
    
    // Convert organelle ID to display name using StringDict
    return StringDict::idToString(organelleId);
}
