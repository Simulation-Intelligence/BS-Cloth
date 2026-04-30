#include "atrpch.h"
#include "Renderer.h"

#include "Loader/Config/Config.h"

namespace ATR
{

    Renderer::Renderer() : 
        config(Config::defaultConfig)
    {  }

    void Renderer::Run()
    {
        try
        {
            InitRenderer();
            InitVulkan();

            while (!this->vkResources.ShouldClose())
                Update();

            Cleanup();
        }
        catch(const Exception& e)
        {
            ATR_ERROR(e.What())
        }
    }

    void Renderer::InitRenderer()
    {
        ATR_LOG_PART("Initializing Renderer");
        ATR_LOG("Running with Config: \n" << this->config);
        this->vkResources.AbsorbConfigs(this->config);
    }

    void Renderer::InitVulkan()
    {
        ATR_LOG_PART("Initializing Vulkan");
        this->vkResources.Init();
    }

    void Renderer::Update()
    {
        this->vkResources.UpdateFrame();
        UpdateStats();
    }

    void Renderer::Cleanup()
    {
        this->vkResources.CleanUp();
    }

    void Renderer::UpdateStats()
    {
        // Disable warning for unused variable. If ATR_VERBOSE is off, the duration variable will be unused
#pragma warning(disable: 4189)      
        auto currentTime = std::chrono::high_resolution_clock::now();
        Float duration = std::chrono::duration<Float, std::chrono::seconds::period>(currentTime - this->lastTime).count();
        this->lastTime = currentTime;
        ATR_LOG_VERBOSE("Rendering Frame " << this->frameCount << " [ FPS: " 
            << std::fixed << std::setfill(' ') << std::right << std::setw(7) << std::setprecision(2) << (1.f / duration) << " ] "
            << this->vkResources.GetUpdateInfo())
        ++frameCount;
#pragma warning(default: 4189)
    }
}
