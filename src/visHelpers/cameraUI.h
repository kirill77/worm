#pragma once

#include "visualization/Window.h"

struct GPUCamera;

class CameraUI
{
    void notifyNewUIState(const UIState& uiState);

private:
    std::shared_ptr<GPUCamera> m_pCamera;
    UIState m_prevUIState;
};

