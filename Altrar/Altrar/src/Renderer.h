#pragma once
#include "atrfwd.h"

#include <chrono>

#include "Loader/Config/Config.h"
#include "VkResources.h"

namespace ATR
{
    class Renderer
    {
    public:
        Renderer();

        void Run();

        // Major Components
        void InitRenderer();
        void InitVulkan();				// Vulkan helpers are all in class `VkResources`
        void Update();
        void Cleanup();

        inline Bool ShouldClose() const { return this->vkResources.ShouldClose(); }

        // Helpers
        void UpdateStats();

        // Proxy: modify mesh
        inline void AddTriangle(std::array<Vertex, 3> vertices) { this->vkResources.AddTriangle(vertices); }
        inline void UpdateMesh(const Mesh& mesh) { this->vkResources.UpdateMesh(mesh); }

        inline void SetScreenshot(Bool ss) { this->config.saveScreenShots = ss; }
        inline void ResumeAtFrame(UInt frame) { this->config.resumeAtFrame = frame; }

    private:
        Config config;
        VkResourceManager vkResources;

        // Performance stats
        UInt frameCount = 0;
        std::chrono::steady_clock::time_point lastTime = std::chrono::high_resolution_clock::now();
    };
}
