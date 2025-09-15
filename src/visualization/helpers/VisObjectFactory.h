#pragma once

#include <memory>
#include "chemistry/molecules/StringDict.h"

// Forward declarations
struct IVisObject;
class Organelle;
struct GPUQueue;

class VisObjectFactory
{
public:
    /**
     * Create appropriate visualization object for the given organelle
     * 
     * @param pOrganelle The organelle to create visualization for
     * @param organelleId The type ID of the organelle
     * @param pQueue GPU queue for resource creation
     * @return Shared pointer to the created visualization object
     */
    static std::shared_ptr<IVisObject> createForOrganelle(
        std::shared_ptr<Organelle> pOrganelle,
        StringDict::ID organelleId,
        GPUQueue* pQueue);
}; 