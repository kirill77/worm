#pragma once

#include <memory>
#include "chemistry/StringDict.h"

// Forward declarations
struct IObjectVis;
class GPUMesh;
class Organelle;
struct Organism;
struct Window;

struct VisObjectContext
{
public:
    std::shared_ptr<IObjectVis> m_pObject;   // the visualized object
    std::shared_ptr<GPUMesh> m_pGpuMesh;     // that object converted to GPUMesh

    // Static factory method to create visualization contexts for organelles
    static std::shared_ptr<VisObjectContext> createForOrganelle(
        std::shared_ptr<Organelle> pOrganelle,
        StringDict::ID organelleId,
        std::shared_ptr<Window> pWindow);
}; 