#pragma once
#include "atrpch.h"

#include "Loader/Config/Config.h"

#include "Geometry/Geometry.h"
#include "VkInfos/VkInfos.h"

namespace ATR
{
    class VkResourceManager
    {
    public:
        VkResourceManager() = default;                  // No explicit constructor needed, init members in declaration

        // Configs
        void AbsorbConfigs(const Config& config);

        // Major Components
        void Init();
        void UpdateFrame();
        inline Bool ShouldClose() const { return glfwWindowShouldClose(this->window); }
        void CleanUp();

        // Initializing Vulkan
        void InitParams();
        void CreateWindow();
        void CreateInstance();
        void SetupDebugMessenger();
        void CreateSurface();
        void SelectPhysicalDevice();
        void CreateLogicalDevice();

        void CreateSwapchain();
        void CreateImageViews();
        void CreateRenderPass();
        void CreateDescriptorSetLayout();
        void CreateGraphicsPipeline();
        void CreateCommandPool();

        void CreateDepthBuffer();
        void CreateTextureImage();
        void CreateVertexBuffer();
        void CreateIndexBuffer();
        void CreateUniformBuffer();
        void CreateFrameBuffers();
        void CreateDescriptorPool();
        void CreateDescriptorSets();
        void CreateCommandBuffer();
        void CreateSyncGadgets();

        // Update
        void DrawFrame();
        void RecreateSwapchain();
        void UpdateImageBuffers();

        // Clean Up
        void CleanUpSwapchain();

        /// Helpers
        // Init
        static VKAPI_ATTR VkBool32 VKAPI_CALL DebugCallback(
            VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity,
            VkDebugUtilsMessageTypeFlagsEXT messageType,
            const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData,
            void* pUserData
        );
        VkResult CreateDebugUtilsMessengerEXT(
            VkInstance instance,
            const VkDebugUtilsMessengerCreateInfoEXT* pCreateInfo,
            const VkAllocationCallbacks* pAllocator,
            VkDebugUtilsMessengerEXT* pDebugMessenger
        );
        VkResult DestroyDebugUtilsMessengerEXT(
            VkInstance instance,
            VkDebugUtilsMessengerEXT debugMessenger,
            const VkAllocationCallbacks* pAllocator
        );

        void GetRequiredExtensions();
        void FindValidationLayers();                            // Validation Layers are specified by the user, not infrastructure
        void RetrieveSwapChainImages();

        Bool DeviceSuitable(VkPhysicalDevice device);
        void FindQueueFamilies(VkPhysicalDevice device);
        Bool CheckDeviceExtensionSupport(VkPhysicalDevice device);
        void QuerySwapChainSupport(VkPhysicalDevice device);
        void ConfigureSwapChain(SwapChainSupportDetails support);
        void CreateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer& buffer, VkDeviceMemory& bufferMemory);
        void CreateStagingBuffer(VkDeviceSize size, VkBuffer& stagingBuffer, VkDeviceMemory& stagingBufferMemory);
        void CopyBuffer(VkBuffer src, VkBuffer dst, VkDeviceSize size);

        void CreateImage(UInt width, UInt height, VkFormat format, VkImageTiling tiling, VkImageUsageFlags usage, VkMemoryPropertyFlags properties, VkImage& image, VkDeviceMemory& imageMemory);
        VkImageView CreateImageView(VkImage image, VkFormat format, VkImageAspectFlags aspectFlags);
        VkFormat FindSupportedFormat(const std::vector<VkFormat>& candidates, VkImageTiling tiling, VkFormatFeatureFlags features);
        inline VkFormat FindDepthFormat();
        inline VkCommandPool& GetTransferCommandPool();
        inline bool HasStencilComponent(VkFormat format);

        VkCommandBuffer CreateCommandBuffer(VkCommandBufferLevel level, VkCommandPool commandPool, Bool startRecording);
        void FlushCommandBuffer(VkCommandBuffer commandBuffer, VkCommandPool commandPool, VkQueue queue, Bool free);
        void InsertImageMemoryBarrier(
            VkCommandBuffer commandBuffer, 
            VkImage image, 
            VkAccessFlags srcAccessMask, 
            VkAccessFlags dstAccessMask, 
            VkImageLayout oldLayout, 
            VkImageLayout newLayout, 
            VkPipelineStageFlags srcStageMask, 
            VkPipelineStageFlags dstStageMask, 
            VkImageSubresourceRange subresourceRange
        );

        void CompileShaders();
        std::vector<char> ReadShaderCode(const char* path);
        VkShaderModule CreateShaderModule(const std::vector<char>& code);

        void SaveScreenshot(String filename, UInt imageIndex);

        // Update
        void RecordCommandBuffer(VkCommandBuffer commandBuffer, UInt imageIndex);
        UInt FindMemoryType(UInt typeFilter, VkMemoryPropertyFlags properties);
        void UpdateUniformBuffer(UInt imageIndex);

        // Getter/Setters
        inline String GetUpdateInfo() { return this->updateInfo; }

        // Proxy
        inline void AddTriangle(std::array<Vertex, 3> vertices) { this->mesh.AddTriangle(vertices); this->meshStale = true; }
        inline void UpdateMesh(const Mesh& mesh) { this->mesh.UpdateMesh(mesh); this->meshStale = true; }

    private:
        // Configs
        UInt width, height;
        Bool enabledValidation;
        String relLocation;
        Bool saveScreenShots;
        std::vector<const char*> instanceExtensions;
        std::vector<const char*> validationLayers;

        const std::vector<const char*> deviceExtensions = {
            VK_KHR_SWAPCHAIN_EXTENSION_NAME
        };

        /// Vulkan Resources
        // Top-level Vulkan Resources
        VkInstance instance;
        VkDebugUtilsMessengerEXT debugMessenger;

        // Vulkan Components
        VkPhysicalDevice physicalDevice;
        VkDevice device;

        VkSurfaceKHR surface;
        VkSwapchainKHR swapchain;
        std::array<VkQueue, QueueFamilyIndices::COUNT> queues;
        std::vector<VkImage> swapchainImages;
        std::vector<VkImageView> swapchainImageViews;
        std::vector<VkFramebuffer> swapchainFrameBuffers;

        VkRenderPass renderPass;
        VkDescriptorSetLayout descriptorSetLayout;
        VkPipelineLayout pipelineLayout;                                // Specify Uniforms
        VkPipeline graphicsPipeline;

        VkCommandPool graphicsCommandPool, transferCommandPool;
        std::vector<VkCommandBuffer> graphicsCommandBuffers;

        VkBuffer vertexBuffer, indexBuffer;                             // Actual buffer on device
        VkDeviceMemory vertexBufferMemory, indexBufferMemory;

        VkImage depthImage;
        VkDeviceMemory depthImageMemory;
        VkImageView depthImageView;

        std::vector<VkBuffer> uniformBuffers;
        std::vector<VkDeviceMemory> uniformBuffersMemory;
        std::vector<void*> uniformBufferMappedMemory;

        VkDescriptorPool descriptorPool;
        std::vector<VkDescriptorSet> descriptorSets;                    // Allocated from descriptorPool, similar to commandBuffers from commandPool

        // Synchronization gadgets
        std::vector<VkSemaphore> imageAvailableSemaphores;
        std::vector<VkSemaphore> renderFinishedSemaphores;
        std::vector<VkFence> inFlightFences;

        // Customized Infos
        QueueFamilyIndices queueIndices;
        SwapChainSupportDetails swapChainSupport;
        SwapChainConfig swapChainConfig;
        
        /// -----------------

        // GLFW Handle
        GLFWwindow* window;

        // Per-update Invariances
        UInt currentFrameNumber = 0;
        UInt currentFrameIndex = 0;
        Bool frameBufferResized = false;
        String updateInfo = "";

        // Static (global) functions
        // TODO may want to have a separate class & files to handle callbacks
        static void FramebufferResizeCallback(GLFWwindow* window, int width, int height);

        // Default Vulkan Configs
        static inline VkDebugUtilsMessengerCreateInfoEXT defaultDebugMessengerCreateInfo = {
            .sType = VK_STRUCTURE_TYPE_DEBUG_UTILS_MESSENGER_CREATE_INFO_EXT,
            .messageSeverity = VK_DEBUG_UTILS_MESSAGE_SEVERITY_VERBOSE_BIT_EXT |
                VK_DEBUG_UTILS_MESSAGE_SEVERITY_WARNING_BIT_EXT |
                VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT,
            .messageType = VK_DEBUG_UTILS_MESSAGE_TYPE_GENERAL_BIT_EXT |
                VK_DEBUG_UTILS_MESSAGE_TYPE_VALIDATION_BIT_EXT |
                VK_DEBUG_UTILS_MESSAGE_TYPE_PERFORMANCE_BIT_EXT,
            .pfnUserCallback = VkResourceManager::DebugCallback,
            .pUserData = nullptr
        };
        static inline constexpr Float defaultQueuePriority = 0.2f;
        static inline constexpr VkClearValue defaultClearValue = { 0.0f, 0.0f, 0.0f, 1.0f };
        static inline constexpr VkClearValue defaultDepthClearValue = {1.f, 0.f};
        static inline constexpr UInt maxFramesInFlight = 2;
        static inline constexpr VkFormat defaultColorFormat = VK_FORMAT_R8G8B8A8_SRGB;

        // Temporary Global Variables
        /*static inline const std::vector<Vertex> vertices = {
            {{-0.5f, -0.5f, 0.f}, {1.0f, 0.0f, 0.0f}},
            {{0.5f, -0.5f, 0.f}, {0.0f, 1.0f, 0.0f}},
            {{0.5f, 0.5f, 0.f}, {0.0f, 0.0f, 1.0f}},
            {{-0.5f, 0.5f, 0.f}, {1.0f, 1.0f, 1.0f}},

            {{-0.5f, -0.5f, -0.5f}, {1.0f, 0.0f, 0.0f}},
            {{0.5f, -0.5f, -0.5f}, {0.0f, 1.0f, 0.0f}},
            {{0.5f, 0.5f, -0.5f}, {0.0f, 0.0f, 1.0f}},
            {{-0.5f, 0.5f, -0.5f}, {1.0f, 1.0f, 1.0f}}
        };

        static inline const std::vector<UInt> indices = {
            0, 1, 2, 2, 3, 0,
            4, 5, 6, 6, 7, 4
        };*/

        Mesh mesh;

        // Update Infos
        bool meshStale = false;
    };

}