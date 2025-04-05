#include "pch.h"
#include "Window.h"
#include "GPUWorld.h"
#include "GPUMesh.h"
#include <iostream>

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow)
{
    // Enable console output for debugging
    AllocConsole();
    FILE* pConsole;
    freopen_s(&pConsole, "CONOUT$", "w", stdout);
    
    std::cout << "Starting visualization application..." << std::endl;
    
    try
    {
        // Create a window
        auto window = std::make_shared<Window>();
        if (!window->createWindowDevicAndSwapChain("TensionSphere Visualization"))
        {
            std::cerr << "Failed to create window and DirectX device." << std::endl;
            return 1;
        }
        
        // Create a GPU world for rendering
        auto world = std::make_shared<GPUWorld>(window);
        
        // Create a simple triangle mesh
        auto mesh = world->createMesh();
        
        // Define vertices
        std::vector<GPUMesh::Vertex> vertices = {
            { float3(0.0f, 0.5f, 0.0f) },   // Top
            { float3(0.5f, -0.5f, 0.0f) },  // Bottom right
            { float3(-0.5f, -0.5f, 0.0f) }, // Bottom left
        };
        
        // Define indices
        std::vector<int3> triangles = {
            int3(0, 1, 2), // Single triangle
        };
        
        // Upload the mesh geometry
        mesh->setGeometry(vertices, triangles);
        
        // Add mesh to the world
        world->addMesh(mesh);
        
        // Main message loop
        MSG msg = {};
        while (msg.message != WM_QUIT)
        {
            // Process window messages
            window->processMessages();
            
            // Render the scene
            world->drawMeshesIntoWindow();
            
            // Handle user input
            const UIState& uiState = window->getCurrentUIState();
            
            // Check for escape key to exit
            if (uiState.getButtonOrKeyPressCount(VK_ESCAPE) > 0)
            {
                PostQuitMessage(0);
            }
            
            // For demo, rotate the camera around the scene
            auto camera = world->getCamera();
            static float angle = 0.0f;
            angle += 0.01f;
            float x = 3.0f * std::sin(angle);
            float z = 3.0f * std::cos(angle);
            camera->setPosition(float3(x, 0.0f, z));
            camera->setLookAt(float3(0.0f, 0.0f, 0.0f));
            
            // Limit frame rate
            Sleep(16); // ~60 FPS
        }
        
        // Clean up
        if (pConsole)
        {
            fclose(pConsole);
        }
        
        return static_cast<int>(msg.wParam);
    }
    catch (const std::exception& e)
    {
        std::cerr << "Exception caught: " << e.what() << std::endl;
        MessageBoxA(nullptr, e.what(), "Error", MB_OK | MB_ICONERROR);
        return 1;
    }
} 