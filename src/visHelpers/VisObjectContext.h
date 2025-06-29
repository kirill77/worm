#pragma once

#include <memory>

// Forward declarations
struct IObjectVis;
class GPUMesh;

struct VisObjectContext
{
    std::shared_ptr<IObjectVis> m_pObject;   // the visualized object
    std::shared_ptr<GPUMesh> m_pGpuMesh;     // that object converted to GPUMesh
}; 