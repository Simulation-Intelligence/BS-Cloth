#include "atrpch.h"
#include "VkResources.h"

#include <chrono>
#include "glm/gtc/matrix_transform.hpp"

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image.h"
#include "stb/stb_image_write.h"

#include "Utility/Profiler.h"

namespace ATR 
{  
    void VkResourceManager::AbsorbConfigs(const Config& config)
    {
        for (const String& layerName : config.validationLayers)
            this->validationLayers.push_back(layerName.c_str());

        this->width = config.width;
        this->height = config.height;
        this->relLocation = config.location;
        this->saveScreenShots = config.saveScreenShots;
        this->currentFrameNumber = config.resumeAtFrame;
    }

    void VkResourceManager::Init()
    {
        this->InitParams();

        // Setup GLFW Window
        this->CreateWindow();

        // Setup Vulkan
        this->CreateInstance();
        this->SetupDebugMessenger();
        this->CreateSurface();
        this->SelectPhysicalDevice();
        this->CreateLogicalDevice();

        // Setup Graphics Pipeline
        this->CreateSwapchain();
        this->CreateImageViews();
        this->CreateRenderPass();
        this->CreateDescriptorSetLayout();
        this->CreateGraphicsPipeline();
        this->CreateCommandPool();

        // Setup Buffers and Syncing
        this->CreateDepthBuffer();
        this->CreateTextureImage();
        this->CreateVertexBuffer();
        this->CreateIndexBuffer();
        this->CreateUniformBuffer();
        this->CreateFrameBuffers();
        this->CreateDescriptorPool();
        this->CreateDescriptorSets();
        this->CreateCommandBuffer();
        this->CreateSyncGadgets();
    }

    void VkResourceManager::UpdateFrame()
    {
        this->updateInfo = "[ Frame Index: " + std::to_string(this->currentFrameIndex) + " ]";
        this->DrawFrame();
        glfwPollEvents();
    }

    void VkResourceManager::CleanUp()
    {
        // Ensure that the device is subject to cleaning up
        vkDeviceWaitIdle(this->device);

        // Clean up validation layer
        if (enabledValidation)
            this->DestroyDebugUtilsMessengerEXT(this->instance, this->debugMessenger, nullptr);

        // Clean up synchronization gadgets
        for (auto& semaphore : this->imageAvailableSemaphores)
            vkDestroySemaphore(this->device, semaphore, nullptr);
        for (auto& semaphore : this->renderFinishedSemaphores)
            vkDestroySemaphore(this->device, semaphore, nullptr);
        for (auto& fence : this->inFlightFences)
            vkDestroyFence(this->device, fence, nullptr);

        // Clean up device-dependent resources
        vkDestroyBuffer(this->device, this->vertexBuffer, nullptr);
        vkFreeMemory(this->device, this->vertexBufferMemory, nullptr);                      // This must be done after the buffer itself has been freed
        vkDestroyBuffer(this->device, this->indexBuffer, nullptr);
        vkFreeMemory(this->device, this->indexBufferMemory, nullptr);

        vkDestroyCommandPool(this->device, this->graphicsCommandPool, nullptr);             // Command buffers are automatically freed when we free the command pool
        vkDestroyCommandPool(this->device, this->transferCommandPool, nullptr);

        vkDestroyPipeline(this->device, this->graphicsPipeline, nullptr);
        vkDestroyRenderPass(this->device, this->renderPass, nullptr);
        vkDestroyPipelineLayout(this->device, this->pipelineLayout, nullptr);
        this->CleanUpSwapchain();

        for (size_t i = 0; i != VkResourceManager::maxFramesInFlight; ++i)
        {
            vkDestroyBuffer(this->device, this->uniformBuffers[i], nullptr);
            vkFreeMemory(this->device, this->uniformBuffersMemory[i], nullptr);
        }

        vkDestroyDescriptorPool(this->device, this->descriptorPool, nullptr);
        vkDestroyDescriptorSetLayout(this->device, this->descriptorSetLayout, nullptr);

        vkDestroyDevice(this->device, nullptr);
        
        // Clean up instance-dependent resources
        vkDestroySurfaceKHR(this->instance, this->surface, nullptr);
        vkDestroyInstance(this->instance, nullptr);

        // Clean up GLFW resources
        glfwDestroyWindow(this->window);
        glfwTerminate();
    }

    void VkResourceManager::InitParams()
    {
        // Prepare for cache folder
        OS::Execute("rmdir /s /q cache");
        OS::Execute("mkdir cache");
    }

    void VkResourceManager::CreateWindow()
    {
        ATR_LOG_PART("Initializing Window");
        glfwInit();
        glfwWindowHint(GLFW_CLIENT_API, GLFW_NO_API);
        glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
        this->window = glfwCreateWindow(this->width, this->height, "Altrar", nullptr, nullptr);
        glfwSetWindowUserPointer(this->window, this);
        glfwSetFramebufferSizeCallback(this->window, VkResourceManager::FramebufferResizeCallback);
        if (!window)
            throw Exception("Failed to create window", ExceptionType::INIT_GLFW);
    }

    void VkResourceManager::CreateInstance()
    {
        ATR_LOG("Creating Vulkan Instance...")
        VkResourceManager::GetRequiredExtensions();
        VkResourceManager::FindValidationLayers();
        
        // Vulkan Instance Creation Info
        // NOTE this cannot be extracted to a function as `appInfo` needs to be available when `VkCreateInstance` is called
        VkApplicationInfo appInfo =
        {
            .sType = VK_STRUCTURE_TYPE_APPLICATION_INFO,
            .pApplicationName = "Altrar",
            .applicationVersion = VK_MAKE_VERSION(1, 0, 0),
            .pEngineName = "No Engine",
            .engineVersion = VK_MAKE_VERSION(1, 0, 0),
            .apiVersion = VK_API_VERSION_1_0
        };

        VkInstanceCreateInfo createInfo =
        {
            .sType = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO,
            .pApplicationInfo = &appInfo,
            .enabledExtensionCount = static_cast<UInt>(this->instanceExtensions.size()),
            .ppEnabledExtensionNames = this->instanceExtensions.data()
        };

        if (enabledValidation)
        {
            createInfo.enabledLayerCount = static_cast<UInt>(this->validationLayers.size());
            createInfo.ppEnabledLayerNames = this->validationLayers.data();

            createInfo.pNext = &VkResourceManager::defaultDebugMessengerCreateInfo;
        }
        else
            createInfo.enabledLayerCount = 0;

        if (vkCreateInstance(&createInfo, nullptr, &this->instance) != VK_SUCCESS)
            throw Exception("Failed to create instance", ExceptionType::INIT_VULKAN);
    }

    void VkResourceManager::SetupDebugMessenger()
    {
        ATR_LOG("Setting Up Debug Structures...")
        if (!this->enabledValidation)
            return;

        if (CreateDebugUtilsMessengerEXT(instance, &VkResourceManager::defaultDebugMessengerCreateInfo, nullptr, &debugMessenger) != VK_SUCCESS)
            throw Exception("Failed to set up debug messenger", ExceptionType::INIT_VULKAN);
    }

    void VkResourceManager::CreateSurface()
    {
        ATR_LOG("Creating Window Surface...")
        if (glfwCreateWindowSurface(this->instance, this->window, nullptr, &this->surface) != VK_SUCCESS)
            throw Exception("Failed to create window surface", ExceptionType::INIT_VULKAN);
    }

    // MSVC Buggy C4702 Detection (unreachable code)
#pragma warning(disable: 4702)
    void VkResourceManager::SelectPhysicalDevice()
    {
        ATR_LOG("Looking for GPUs...")
        this->physicalDevice = VK_NULL_HANDLE;
        UInt physicalDeviceCount = 0;
        vkEnumeratePhysicalDevices(this->instance, &physicalDeviceCount, nullptr);

        std::vector<VkPhysicalDevice> availablePhysicalDevices(physicalDeviceCount);
        vkEnumeratePhysicalDevices(this->instance, &physicalDeviceCount, availablePhysicalDevices.data());

        if (availablePhysicalDevices.size() == 0)
            throw Exception("No available GPU Found.", ExceptionType::INIT_VULKAN);

        for (auto& device : availablePhysicalDevices)
        {
            if (DeviceSuitable(device))
            {
                this->physicalDevice = device;
                break;
            }
            throw Exception("No Suitable GPU Found.", ExceptionType::INIT_VULKAN);
        }

        VkPhysicalDeviceProperties deviceProperties;
        vkGetPhysicalDeviceProperties(this->physicalDevice, &deviceProperties);
        ATR_PRINT_VERBOSE("Using Physical Device: " + String(deviceProperties.deviceName))
        ATR_PRINT_VERBOSE("Queue Family Indices: \n" << this->queueIndices)
    }
#pragma warning(default: 4702)

    void VkResourceManager::CreateLogicalDevice()
    {
        ATR_LOG("Setting Up Device...")
        VkPhysicalDeviceFeatures deviceFeatures = {};

        std::vector<VkDeviceQueueCreateInfo> queueCreateInfos;
        std::set<UInt> uniqueQueueFamilies;
        for (size_t index = 0; index != QueueFamilyIndices::COUNT; ++index)
            uniqueQueueFamilies.insert(this->queueIndices.indices[index].value());

        for (UInt uniqueQueueIndex : uniqueQueueFamilies)
        {
            VkDeviceQueueCreateInfo queueCreateInfo = {
                .sType = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO,
                //.queueFamilyIndex = this->queueIndices.indices[index].value(),
                .queueFamilyIndex = uniqueQueueIndex,
                .queueCount = 1,
                .pQueuePriorities = &VkResourceManager::defaultQueuePriority
            };
            queueCreateInfos.push_back(queueCreateInfo);
        }

        VkDeviceCreateInfo deviceCreateInfo = {
            .sType = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO,
            .queueCreateInfoCount = static_cast<UInt>(queueCreateInfos.size()),
            .pQueueCreateInfos = queueCreateInfos.data(),
            .enabledExtensionCount = (UInt)(this->deviceExtensions.size()),
            .ppEnabledExtensionNames = this->deviceExtensions.data(),
            .pEnabledFeatures = &deviceFeatures
        };

        // NOTE from Vulkan 1.3.290 this is not necessary as device will automatically have the same validation layers as the instance
        if (this->enabledValidation)
        {
            deviceCreateInfo.enabledLayerCount = static_cast<UInt>(this->validationLayers.size());
            deviceCreateInfo.ppEnabledLayerNames = this->validationLayers.data();
        }
        else
            deviceCreateInfo.enabledLayerCount = 0;

        VkResult result = vkCreateDevice(this->physicalDevice, &deviceCreateInfo, nullptr, &this->device);
        if (result != VK_SUCCESS)
            throw Exception("Failed to create logical device", ExceptionType::INIT_VULKAN);

        for (size_t index = 0; index != QueueFamilyIndices::COUNT; ++index)
            vkGetDeviceQueue(this->device, this->queueIndices.indices[index].value(), 0, &this->queues[index]);
    }

    void VkResourceManager::CreateSwapchain()
    {
        ATR_LOG("Creating Swapchain...")
        this->QuerySwapChainSupport(this->physicalDevice);
        this->ConfigureSwapChain(this->swapChainSupport);

        UInt imageCount = this->swapChainSupport.capabilities.minImageCount + 1;
        UInt maxSuportedImages = this->swapChainSupport.capabilities.maxImageCount;
        if (maxSuportedImages != 0 && imageCount > maxSuportedImages)
            imageCount = maxSuportedImages;

        VkSwapchainCreateInfoKHR createInfo = {
            .sType = VK_STRUCTURE_TYPE_SWAPCHAIN_CREATE_INFO_KHR,
            .surface = this->surface,
            .minImageCount = imageCount,
            .imageFormat = this->swapChainConfig.format.format,
            .imageColorSpace = this->swapChainConfig.format.colorSpace,
            .imageExtent = this->swapChainConfig.extent,
            .imageArrayLayers = 1,
            .imageUsage = VK_IMAGE_USAGE_COLOR_ATTACHMENT_BIT | VK_IMAGE_USAGE_TRANSFER_SRC_BIT,
            .preTransform = this->swapChainSupport.capabilities.currentTransform,
            .compositeAlpha = VK_COMPOSITE_ALPHA_OPAQUE_BIT_KHR,
            .presentMode = this->swapChainConfig.presentMode,
            .clipped = VK_TRUE,
            .oldSwapchain = VK_NULL_HANDLE
        };

        std::vector<UInt> uniqueIndices = this->queueIndices.UniqueIndices();
        if (!this->queueIndices.AllQueuesSame())
        {
            createInfo.imageSharingMode = VK_SHARING_MODE_CONCURRENT;
            createInfo.queueFamilyIndexCount = static_cast<UInt>(uniqueIndices.size());
            createInfo.pQueueFamilyIndices = uniqueIndices.data();
        }
        else
            createInfo.imageSharingMode = VK_SHARING_MODE_EXCLUSIVE;

        if (vkCreateSwapchainKHR(this->device, &createInfo, nullptr, &this->swapchain) != VK_SUCCESS)
            throw Exception("Failed to create swapchain", ExceptionType::INIT_VULKAN);

        this->RetrieveSwapChainImages();
    }

    void VkResourceManager::CreateImageViews()
    {
        ATR_LOG("Creating image views from swapchain...")
        swapchainImageViews.resize(swapchainImages.size());
        for (size_t i = 0; i != swapchainImages.size(); ++i)
        {
            // TODO replace this with CreateImageView
            VkImageViewCreateInfo createInfo = {
                .sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO,
                .image = swapchainImages[i],
                .viewType = VK_IMAGE_VIEW_TYPE_2D,
                .format = swapChainConfig.format.format,

                .components = {
                    .r = VK_COMPONENT_SWIZZLE_IDENTITY,
                    .g = VK_COMPONENT_SWIZZLE_IDENTITY,
                    .b = VK_COMPONENT_SWIZZLE_IDENTITY,
                    .a = VK_COMPONENT_SWIZZLE_IDENTITY
                },

                .subresourceRange = {
                    .aspectMask = VK_IMAGE_ASPECT_COLOR_BIT,
                    .baseMipLevel = 0,
                    .levelCount = 1,
                    .baseArrayLayer = 0,
                    .layerCount = 1
                }
            };

            if (vkCreateImageView(this->device, &createInfo, nullptr, &this->swapchainImageViews[i]) != VK_SUCCESS)
                throw Exception("Failed to create swapchain image views", ExceptionType::INIT_VULKAN);
        }
    }

    void VkResourceManager::CreateRenderPass()
    {
        ATR_LOG("Creating Render Pass...")

        VkAttachmentDescription colorAttachment = {
            .format = this->swapChainConfig.format.format,
            .samples = VK_SAMPLE_COUNT_1_BIT,
            .loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR,
            .storeOp = VK_ATTACHMENT_STORE_OP_STORE,
            .stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE,
            .stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
            .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
            .finalLayout = VK_IMAGE_LAYOUT_PRESENT_SRC_KHR
        };

        VkAttachmentReference colorAttachmentRef = {
            .attachment = 0,
            .layout = VK_IMAGE_LAYOUT_COLOR_ATTACHMENT_OPTIMAL
        };

        VkAttachmentDescription depthAttachment = {
            .format = this->FindDepthFormat(),
            .samples = VK_SAMPLE_COUNT_1_BIT,
            .loadOp = VK_ATTACHMENT_LOAD_OP_CLEAR,
            .storeOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
            .stencilLoadOp = VK_ATTACHMENT_LOAD_OP_DONT_CARE,
            .stencilStoreOp = VK_ATTACHMENT_STORE_OP_DONT_CARE,
            .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED,
            .finalLayout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL
        };

        VkAttachmentReference depthAttachmentRef = {
            .attachment = 1,
            .layout = VK_IMAGE_LAYOUT_DEPTH_STENCIL_ATTACHMENT_OPTIMAL
        };

        VkSubpassDependency dependency = {
            .srcSubpass = VK_SUBPASS_EXTERNAL,
            .dstSubpass = 0,
            .srcStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT,
            .dstStageMask = VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT | VK_PIPELINE_STAGE_EARLY_FRAGMENT_TESTS_BIT,
            .srcAccessMask = 0,
            .dstAccessMask = VK_ACCESS_COLOR_ATTACHMENT_WRITE_BIT | VK_ACCESS_DEPTH_STENCIL_ATTACHMENT_WRITE_BIT
        };

        VkSubpassDescription subpass = {
            .pipelineBindPoint = VK_PIPELINE_BIND_POINT_GRAPHICS,
            .colorAttachmentCount = 1,
            .pColorAttachments = &colorAttachmentRef,
            .pDepthStencilAttachment = &depthAttachmentRef
        };

        std::array<VkAttachmentDescription, 2> attachments = { colorAttachment, depthAttachment };

        VkRenderPassCreateInfo renderPass = {
            .sType = VK_STRUCTURE_TYPE_RENDER_PASS_CREATE_INFO,
            .attachmentCount = static_cast<UInt>(attachments.size()),
            .pAttachments = attachments.data(),
            .subpassCount = 1,
            .pSubpasses = &subpass,
            .dependencyCount = 1,
            .pDependencies = &dependency
        };

        if (vkCreateRenderPass(this->device, &renderPass, nullptr, &this->renderPass) != VK_SUCCESS)
            throw Exception("Failed to create render pass", ExceptionType::INIT_PIPELINE);
    }

    void VkResourceManager::CreateDescriptorSetLayout()
    {
        ATR_LOG("Creating Descriptor Set Layout...")
        VkDescriptorSetLayoutBinding uboLayoutBinding = {
            .binding = 0,
            .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
            .descriptorCount = 1,
            .stageFlags = VK_SHADER_STAGE_VERTEX_BIT,
            .pImmutableSamplers = nullptr
        };

        VkDescriptorSetLayoutCreateInfo layoutInfo = {
            .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO,
            .bindingCount = 1,
            .pBindings = &uboLayoutBinding
        };

        if (vkCreateDescriptorSetLayout(this->device, &layoutInfo, nullptr, &this->descriptorSetLayout) != VK_SUCCESS)
            throw Exception("Failed to create descriptor set layout", ExceptionType::INIT_PIPELINE);
    }

    void VkResourceManager::CreateGraphicsPipeline()
    {
        ATR_LOG("Creating Graphics Pipeline...")

        // Compile Shaders
        // TODO refactor and define paths in ATRConst.h
        ATR_LOG_SUB("Compiling Shaders...")
        this->CompileShaders();

        std::vector<char> vertShaderCode = ReadShaderCode("bin/shaders/vert.spv");
        std::vector<char> fragShaderCode = ReadShaderCode("bin/shaders/frag.spv");

        // Shader modules can be destroyed after shader stage creation, and therefore is not a member variable of VkResourceManager
        VkShaderModule vertShaderModule = CreateShaderModule(vertShaderCode);
        VkShaderModule fragShaderModule = CreateShaderModule(fragShaderCode);

        // Shader Stage Creation
        ATR_LOG_SUB("Creating Shader Stages...")
        VkPipelineShaderStageCreateInfo vertShaderStageInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
            .stage = VK_SHADER_STAGE_VERTEX_BIT,
            .module = vertShaderModule,
            .pName = Const::DefaultShaderEntryPoint
        };
        VkPipelineShaderStageCreateInfo fragShaderStageInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO,
            .stage = VK_SHADER_STAGE_FRAGMENT_BIT,
            .module = fragShaderModule,
            .pName = Const::DefaultShaderEntryPoint
        };
        VkPipelineShaderStageCreateInfo shaderStages[] = { vertShaderStageInfo, fragShaderStageInfo };

        // Configurable states
        std::vector<VkDynamicState> dynamicStates = {
            VK_DYNAMIC_STATE_VIEWPORT,
            VK_DYNAMIC_STATE_SCISSOR
        };

        VkPipelineDynamicStateCreateInfo dynamicStateCreateInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_DYNAMIC_STATE_CREATE_INFO,
            .dynamicStateCount = static_cast<UInt>(dynamicStates.size()),
            .pDynamicStates = dynamicStates.data()
        };

        // Vertex Layout: TODO move the hardcoded values in the shader to the configurable ones here
        ATR_LOG_SUB("Configuring Input Layouts...")
        VkPipelineVertexInputStateCreateInfo vertexInputInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_VERTEX_INPUT_STATE_CREATE_INFO,
            .vertexBindingDescriptionCount = 1,
            .pVertexBindingDescriptions = &Vertex::bindingDescription,
            .vertexAttributeDescriptionCount = static_cast<UInt>(Vertex::attributeDescriptions.size()),
            .pVertexAttributeDescriptions = Vertex::attributeDescriptions.data(),
        };

        // Input Assembly
        VkPipelineInputAssemblyStateCreateInfo inputAssemblyInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_INPUT_ASSEMBLY_STATE_CREATE_INFO,
            .topology = VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST,
            .primitiveRestartEnable = VK_FALSE
        };

        // Viewport and Scissor (Viewport States)
        ATR_LOG_SUB("Configuring Viewport...")
        VkViewport viewport = {
            .x = 0.0f,
            .y = 0.0f,
            .width = (Float)this->swapChainConfig.extent.width,
            .height = (Float)this->swapChainConfig.extent.height,
            .minDepth = 0.0f,           // Notice Vulkan uses depth from 0 to 1
            .maxDepth = 1.0f
        };

        VkRect2D scissor = {
            .offset = { 0, 0 },
            .extent = this->swapChainConfig.extent
        };

        // TODO see how this coops with dynamic configuration
        VkPipelineViewportStateCreateInfo viewportStateInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_VIEWPORT_STATE_CREATE_INFO,
            .viewportCount = 1,
            .pViewports = &viewport,
            .scissorCount = 1,
            .pScissors = &scissor
        };

        // Rasterizer
        ATR_LOG_SUB("Configuring Rasterizer...")
        VkPipelineRasterizationStateCreateInfo rasterizerInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_RASTERIZATION_STATE_CREATE_INFO,
            .depthClampEnable = VK_FALSE,
            .rasterizerDiscardEnable = VK_FALSE,
            .polygonMode = VK_POLYGON_MODE_FILL,
            .cullMode = VK_CULL_MODE_NONE,
            .frontFace = VK_FRONT_FACE_COUNTER_CLOCKWISE,            
            .depthBiasEnable = VK_FALSE,
            .lineWidth = 1.0f
        };

        // Multisampling
        VkPipelineMultisampleStateCreateInfo multisamplingInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_MULTISAMPLE_STATE_CREATE_INFO,
            .rasterizationSamples = VK_SAMPLE_COUNT_1_BIT,
            .sampleShadingEnable = VK_FALSE,
            .minSampleShading = 1.0f,
            .pSampleMask = nullptr,
            .alphaToCoverageEnable = VK_FALSE,
            .alphaToOneEnable = VK_FALSE
        };

        VkPipelineDepthStencilStateCreateInfo depthStencilInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_DEPTH_STENCIL_STATE_CREATE_INFO,
            .depthTestEnable = VK_TRUE,
            .depthWriteEnable = VK_TRUE,
            .depthCompareOp = VK_COMPARE_OP_LESS,
            .depthBoundsTestEnable = VK_FALSE,
            .stencilTestEnable = VK_FALSE,
            .front = {},
            .back = {},
            .minDepthBounds = 0.0f,
            .maxDepthBounds = 1.0f,
        };

        // Color Blending per buffer
        VkPipelineColorBlendAttachmentState colorBlendAttachment = {
            // Enable Alpha Blending
            .blendEnable = VK_TRUE,
            .srcColorBlendFactor = VK_BLEND_FACTOR_SRC_ALPHA,
            .dstColorBlendFactor = VK_BLEND_FACTOR_ONE_MINUS_SRC_ALPHA,
            .colorBlendOp = VK_BLEND_OP_ADD,
            .srcAlphaBlendFactor = VK_BLEND_FACTOR_ONE,
            .dstAlphaBlendFactor = VK_BLEND_FACTOR_ZERO,
            .alphaBlendOp = VK_BLEND_OP_ADD,
            .colorWriteMask = VK_COLOR_COMPONENT_R_BIT | VK_COLOR_COMPONENT_G_BIT | VK_COLOR_COMPONENT_B_BIT | VK_COLOR_COMPONENT_A_BIT
        };

        // Color Blending per Pipeline
        VkPipelineColorBlendStateCreateInfo colorBlendingInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_COLOR_BLEND_STATE_CREATE_INFO,
            .logicOpEnable = VK_FALSE,
            .logicOp = VK_LOGIC_OP_COPY,
            .attachmentCount = 1,
            .pAttachments = &colorBlendAttachment,
            .blendConstants = { 0.0f, 0.0f, 0.0f, 0.0f }
        };

        // Pipeline Layout for Uniforms
        VkPipelineLayoutCreateInfo pipelineLayoutInfo = {
            .sType = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO,
            .setLayoutCount = 1,
            .pSetLayouts = &this->descriptorSetLayout,
            .pushConstantRangeCount = 0,
            .pPushConstantRanges = nullptr
        };

        if (vkCreatePipelineLayout(this->device, &pipelineLayoutInfo, nullptr, &this->pipelineLayout) != VK_SUCCESS)
            throw Exception("Failed to create pipeline layout", ExceptionType::INIT_PIPELINE);

        // Creating the Graphics Pipeline
        VkGraphicsPipelineCreateInfo createInfo = {
            .sType = VK_STRUCTURE_TYPE_GRAPHICS_PIPELINE_CREATE_INFO,
            .stageCount = 2,
            .pStages = shaderStages,
            .pVertexInputState = &vertexInputInfo,
            .pInputAssemblyState = &inputAssemblyInfo,
            .pViewportState = &viewportStateInfo,
            .pRasterizationState = &rasterizerInfo,
            .pMultisampleState = &multisamplingInfo,
            .pDepthStencilState = &depthStencilInfo,
            .pColorBlendState = &colorBlendingInfo,
            .pDynamicState = &dynamicStateCreateInfo,
            .layout = this->pipelineLayout,
            .renderPass = this->renderPass, 
            .subpass = 0,
            .basePipelineHandle = VK_NULL_HANDLE,
            .basePipelineIndex = -1
        };

        if (vkCreateGraphicsPipelines(this->device, VK_NULL_HANDLE, 1, &createInfo, nullptr, &this->graphicsPipeline) != VK_SUCCESS)
            throw Exception("Failed to create graphics pipeline", ExceptionType::INIT_PIPELINE);

        ATR_LOG("Graphics Pipeline Created Successfully.")

        // Cleanup Loaded ShaderCode
        vkDestroyShaderModule(this->device, vertShaderModule, nullptr);
        vkDestroyShaderModule(this->device, fragShaderModule, nullptr);
    }

    void VkResourceManager::CreateFrameBuffers()
    {
        ATR_LOG("Creating Framebuffers...")

        this->swapchainFrameBuffers.resize(this->swapchainImageViews.size());

        for (size_t i = 0; i != this->swapchainImageViews.size(); ++i)
        {
            std::array<VkImageView, 2> attachments = { this->swapchainImageViews[i], this->depthImageView };

            VkFramebufferCreateInfo framebufferInfo = {
                .sType = VK_STRUCTURE_TYPE_FRAMEBUFFER_CREATE_INFO,
                .renderPass = this->renderPass,
                .attachmentCount = static_cast<UInt>(attachments.size()),
                .pAttachments = attachments.data(),
                .width = this->swapChainConfig.extent.width,
                .height = this->swapChainConfig.extent.height,
                .layers = 1
            };

            if (vkCreateFramebuffer(this->device, &framebufferInfo, nullptr, &this->swapchainFrameBuffers[i]) != VK_SUCCESS)
                throw Exception("Failed to create framebuffer", ExceptionType::INIT_PIPELINE);
        }
    }

    void VkResourceManager::CreateCommandPool()
    {
        ATR_LOG("Creating Command Pool...")
        VkCommandPoolCreateInfo graphicsCommandPoolInfo = {
            .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
            .flags = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT,
            .queueFamilyIndex = this->queueIndices.indices[QueueFamilyIndices::GRAPHICS].value()
        };

        if (vkCreateCommandPool(this->device, &graphicsCommandPoolInfo, nullptr, &this->graphicsCommandPool) != VK_SUCCESS)
            throw Exception("Failed to create graphics command pool", ExceptionType::INIT_PIPELINE);

        // if the graphics and transfer queue families are different, we also need a separate pool for transfer commands
        if (this->queueIndices.SeparateTransferQueue())
        {
            VkCommandPoolCreateInfo transferCommandPoolInfo = {
                .sType = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO,
                .flags = VK_COMMAND_POOL_CREATE_TRANSIENT_BIT,
                .queueFamilyIndex = this->queueIndices.indices[QueueFamilyIndices::TRANSFER].value()
            };

            if (vkCreateCommandPool(this->device, &transferCommandPoolInfo, nullptr, &this->transferCommandPool) != VK_SUCCESS)
                throw Exception("Failed to create transfer command pool", ExceptionType::INIT_PIPELINE);
        }
    }

    void VkResourceManager::CreateDepthBuffer()
    {
        VkFormat depthFormat = this->FindDepthFormat();
        this->CreateImage(
            this->swapChainConfig.extent.width,
            this->swapChainConfig.extent.height,
            depthFormat,
            VK_IMAGE_TILING_OPTIMAL,
            VK_IMAGE_USAGE_DEPTH_STENCIL_ATTACHMENT_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            this->depthImage,
            this->depthImageMemory
        );
        this->depthImageView = this->CreateImageView(this->depthImage, depthFormat, VK_IMAGE_ASPECT_DEPTH_BIT);
    }

    void VkResourceManager::CreateTextureImage()
    {
        // TODO
    }

    void VkResourceManager::CreateVertexBuffer()
    {
        // TODO make this more flexible
        ATR_LOG("Creating Vertex Buffer...")
        VkDeviceSize bufferSize = sizeof(this->mesh.GetVertices()[0]) * this->mesh.GetVertices().size();

        VkBuffer stagingBuffer;                         // Buffer on CPU, temporary, host-visible
        VkDeviceMemory stagingBufferMemory;

        this->CreateStagingBuffer(bufferSize, stagingBuffer, stagingBufferMemory);

        void* data;
        vkMapMemory(this->device, stagingBufferMemory, 0, bufferSize, 0, &data);
            memcpy(data, this->mesh.GetVertices().data(), static_cast<size_t>(bufferSize));
        vkUnmapMemory(this->device, stagingBufferMemory);

        this->CreateBuffer(
            bufferSize,
            VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,       // Will receive transfer from the host, and used as vertex buffer
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,                                        // Most suitable (performative) for device local access
            this->vertexBuffer,
            this->vertexBufferMemory
        );

        this->CopyBuffer(stagingBuffer, this->vertexBuffer, bufferSize);

        vkDestroyBuffer(this->device, stagingBuffer, nullptr);
        vkFreeMemory(this->device, stagingBufferMemory, nullptr);
    }

    void VkResourceManager::CreateIndexBuffer()
    {
        ATR_LOG("Creating Index Buffer...")
        VkDeviceSize bufferSize = sizeof(this->mesh.GetIndices()[0]) * this->mesh.GetIndices().size();

        VkBuffer stagingBuffer;
        VkDeviceMemory stagingBufferMemory;

        this->CreateStagingBuffer(bufferSize, stagingBuffer, stagingBufferMemory);

        void* data;
        vkMapMemory(this->device, stagingBufferMemory, 0, bufferSize, 0, &data);
            memcpy(data, this->mesh.GetIndices().data(), static_cast<size_t>(bufferSize));
        vkUnmapMemory(this->device, stagingBufferMemory);

        this->CreateBuffer(
            bufferSize,
            VK_BUFFER_USAGE_TRANSFER_DST_BIT | VK_BUFFER_USAGE_INDEX_BUFFER_BIT,
            VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
            this->indexBuffer,
            this->indexBufferMemory
        );

        this->CopyBuffer(stagingBuffer, this->indexBuffer, bufferSize);

        vkDestroyBuffer(this->device, stagingBuffer, nullptr);
        vkFreeMemory(this->device, stagingBufferMemory, nullptr);
    }

    void VkResourceManager::CreateUniformBuffer()
    {
        ATR_LOG("Creating Uniform Buffer...")
        VkDeviceSize bufferSize = sizeof(UniformBufferObject);

        this->uniformBuffers.resize(VkResourceManager::maxFramesInFlight);
        this->uniformBuffersMemory.resize(VkResourceManager::maxFramesInFlight);
        this->uniformBufferMappedMemory.resize(VkResourceManager::maxFramesInFlight);

        // We are not using staging buffer here because the uniform buffer is updated every frame
        //  Frequent allocation may in fact hamper performance
        for (UInt i = 0; i != VkResourceManager::maxFramesInFlight; ++i)
        {
            this->CreateBuffer(
                bufferSize,
                VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT,
                VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                this->uniformBuffers[i],
                this->uniformBuffersMemory[i]
            );

            // no vkUnmap is here as uniforms may be updated frequently throughout the application
            //  a persistent memory mapping is necessary (and more performant)
            vkMapMemory(this->device, this->uniformBuffersMemory[i], 0, bufferSize, 0, &this->uniformBufferMappedMemory[i]);
        }
    }

    void VkResourceManager::CreateDescriptorPool()
    {
        ATR_LOG("Creating Descriptor Pool...")
        VkDescriptorPoolSize poolSize = {
            .type = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
            .descriptorCount = static_cast<UInt>(VkResourceManager::maxFramesInFlight)
        };

        VkDescriptorPoolCreateInfo poolInfo = {
            .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO,
            .maxSets = static_cast<UInt>(VkResourceManager::maxFramesInFlight),
            .poolSizeCount = 1,
            .pPoolSizes = &poolSize
        };

        if (vkCreateDescriptorPool(this->device, &poolInfo, nullptr, &this->descriptorPool) != VK_SUCCESS)
            throw Exception("Failed to create descriptor pool", ExceptionType::INIT_BUFFER);
    }

    void VkResourceManager::CreateDescriptorSets()
    {
        ATR_LOG("Creating Descriptor Sets...")
        std::vector<VkDescriptorSetLayout> layouts(VkResourceManager::maxFramesInFlight, this->descriptorSetLayout);

        // One descriptor set for each frame in flight, with the same layout (the one specified by struct UniformBufferObject)
        VkDescriptorSetAllocateInfo allocInfo = {
            .sType = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO,
            .descriptorPool = this->descriptorPool,
            .descriptorSetCount = static_cast<UInt>(VkResourceManager::maxFramesInFlight),
            .pSetLayouts = layouts.data()
        };

        this->descriptorSets.resize(VkResourceManager::maxFramesInFlight);
        if (vkAllocateDescriptorSets(this->device, &allocInfo, this->descriptorSets.data()) != VK_SUCCESS)
            throw Exception("Failed to allocate descriptor sets", ExceptionType::INIT_BUFFER);

        // Now fill in the descriptor sets
        for (UInt i = 0; i != VkResourceManager::maxFramesInFlight; ++i)
        {
            VkDescriptorBufferInfo bufferInfo = {
                .buffer = this->uniformBuffers[i],
                .offset = 0,
                .range = sizeof(UniformBufferObject)
            };

            VkWriteDescriptorSet descriptorWrite = {
                .sType = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET,
                .dstSet = this->descriptorSets[i],
                .dstBinding = 0,
                .dstArrayElement = 0,
                .descriptorCount = 1,
                .descriptorType = VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER,
                .pBufferInfo = &bufferInfo
            };

            vkUpdateDescriptorSets(this->device, 1, &descriptorWrite, 0, nullptr);
        }
    }

    void VkResourceManager::CreateCommandBuffer()
    {
        ATR_LOG("Creating Command Buffer...")
        this->graphicsCommandBuffers.resize(VkResourceManager::maxFramesInFlight);
        VkCommandBufferAllocateInfo commandBufferInfo = {
            .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
            .commandPool = this->graphicsCommandPool,
            .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
            .commandBufferCount = static_cast<UInt>(this->graphicsCommandBuffers.size())
        };

        if (vkAllocateCommandBuffers(this->device, &commandBufferInfo, graphicsCommandBuffers.data()) != VK_SUCCESS)
            throw Exception("Failed to allocate command buffers", ExceptionType::INIT_BUFFER);
    }

    void VkResourceManager::CreateSyncGadgets()
    {
        ATR_LOG("Creating Synchronization Gadgets...")

        this->imageAvailableSemaphores.resize(VkResourceManager::maxFramesInFlight);
        this->renderFinishedSemaphores.resize(VkResourceManager::maxFramesInFlight);
        this->inFlightFences.resize(VkResourceManager::maxFramesInFlight);

        VkSemaphoreCreateInfo semaphoreInfo = {
            .sType = VK_STRUCTURE_TYPE_SEMAPHORE_CREATE_INFO
        };
        VkFenceCreateInfo fenceInfo = {
            .sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO,
            .flags = VK_FENCE_CREATE_SIGNALED_BIT
        };

        for (UInt i = 0; i != VkResourceManager::maxFramesInFlight; ++i)
        {
            if (vkCreateSemaphore(this->device, &semaphoreInfo, nullptr, &this->imageAvailableSemaphores[i]) != VK_SUCCESS ||
                vkCreateSemaphore(this->device, &semaphoreInfo, nullptr, &this->renderFinishedSemaphores[i]) != VK_SUCCESS)
                throw Exception("Failed to create synchronization gadgets: semaphores", ExceptionType::INIT_PIPELINE);

            if (vkCreateFence(this->device, &fenceInfo, nullptr, &this->inFlightFences[i]) != VK_SUCCESS)
                throw Exception("Failed to create synchronization gadgets: fence", ExceptionType::INIT_PIPELINE);
        }
    }

    void VkResourceManager::DrawFrame()
    {
        if (this->meshStale)
        {
            this->UpdateImageBuffers();
            this->meshStale = false;
        }

        vkWaitForFences(this->device, 1, &this->inFlightFences[this->currentFrameIndex], VK_TRUE, UINT64_MAX);

        UInt imageIndex;
        VkResult result = vkAcquireNextImageKHR(this->device, this->swapchain, UINT64_MAX, this->imageAvailableSemaphores[this->currentFrameIndex], VK_NULL_HANDLE, &imageIndex);

        if (result == VK_ERROR_OUT_OF_DATE_KHR)
        {
            this->RecreateSwapchain();
            return;
        }
        else if (result != VK_SUCCESS && result != VK_SUBOPTIMAL_KHR)
            throw Exception("Failed to acquire swapchain image", ExceptionType::UPDATE_RENDER);

        vkResetFences(this->device, 1, &this->inFlightFences[this->currentFrameIndex]);             // Reset here to avoid deadlock

        vkResetCommandBuffer(this->graphicsCommandBuffers[this->currentFrameIndex], 0);
        RecordCommandBuffer(this->graphicsCommandBuffers[this->currentFrameIndex], imageIndex);

        VkSemaphore waitSemaphores[] = { this->imageAvailableSemaphores[this->currentFrameIndex] };
        VkSemaphore signalSemaphores[] = { this->renderFinishedSemaphores[this->currentFrameIndex] };
        VkPipelineStageFlags waitStages[] = { VK_PIPELINE_STAGE_COLOR_ATTACHMENT_OUTPUT_BIT };

        this->UpdateUniformBuffer(this->currentFrameIndex);

        VkSubmitInfo submitInfo = {
            .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
            .waitSemaphoreCount = 1,
            .pWaitSemaphores = waitSemaphores,
            .pWaitDstStageMask = waitStages,
            .commandBufferCount = 1,
            .pCommandBuffers = &this->graphicsCommandBuffers[this->currentFrameIndex],
            .signalSemaphoreCount = 1,
            .pSignalSemaphores = &this->renderFinishedSemaphores[this->currentFrameIndex]
        };

        if (vkQueueSubmit(this->queues[QueueFamilyIndices::GRAPHICS], 1, &submitInfo, this->inFlightFences[this->currentFrameIndex]) != VK_SUCCESS)
            throw Exception("Failed to submit draw command buffer", ExceptionType::UPDATE_RENDER);

        VkSwapchainKHR swapchains[] = { this->swapchain };
        VkPresentInfoKHR presentInfo = {
            .sType = VK_STRUCTURE_TYPE_PRESENT_INFO_KHR,
            .waitSemaphoreCount = 1,
            .pWaitSemaphores = signalSemaphores,
            .swapchainCount = 1,
            .pSwapchains = swapchains,
            .pImageIndices = &imageIndex,
            .pResults = nullptr
        };

        result = vkQueuePresentKHR(this->queues[QueueFamilyIndices::PRESENT], &presentInfo);
        if (result == VK_ERROR_OUT_OF_DATE_KHR || result == VK_SUBOPTIMAL_KHR || this->frameBufferResized)
        {
            this->RecreateSwapchain();
            this->frameBufferResized = false;
        }
        else if (result != VK_SUCCESS)
            throw Exception("Failed to present swapchain image", ExceptionType::UPDATE_RENDER);

        if (this->saveScreenShots)
            this->SaveScreenshot(String("cache/frame_") + std::to_string(this->currentFrameNumber) + ".png", imageIndex);

        ++this->currentFrameNumber;
        this->currentFrameIndex = (this->currentFrameIndex + 1) % VkResourceManager::maxFramesInFlight;
    }

    void VkResourceManager::CleanUpSwapchain()
    {
        for (const auto& view : this->swapchainImageViews)
            vkDestroyImageView(this->device, view, nullptr);
        for (const auto& framebuffer : this->swapchainFrameBuffers)
            vkDestroyFramebuffer(this->device, framebuffer, nullptr);

        vkDestroyImageView(this->device, this->depthImageView, nullptr);
        vkDestroyImage(this->device, this->depthImage, nullptr);
        vkFreeMemory(this->device, this->depthImageMemory, nullptr);

        vkDestroySwapchainKHR(this->device, this->swapchain, nullptr);
    }

    // Function signature must conform to Vulkan, so there can be unused variables
#pragma warning(disable: 4100)
    VKAPI_ATTR VkBool32 VKAPI_CALL VkResourceManager::DebugCallback(\
        VkDebugUtilsMessageSeverityFlagBitsEXT messageSeverity, 
        VkDebugUtilsMessageTypeFlagsEXT messageType, 
        const VkDebugUtilsMessengerCallbackDataEXT* pCallbackData, 
        void* pUserData
    )
    {
        if (messageSeverity >= VK_DEBUG_UTILS_MESSAGE_SEVERITY_ERROR_BIT_EXT)
            ATR_ERROR(String("[Validation Layer] ") + pCallbackData->pMessage)
        else
            ATR_PRINT_VERBOSE(String("[Validation Layer] ") + pCallbackData->pMessage)

        return VK_FALSE;
    }
#pragma warning(default: 4100)

    VkResult VkResourceManager::CreateDebugUtilsMessengerEXT(VkInstance instance, 
        const VkDebugUtilsMessengerCreateInfoEXT* pCreateInfo, 
        const VkAllocationCallbacks* pAllocator, 
        VkDebugUtilsMessengerEXT* pDebugMessenger
    )
    {
        auto func = (PFN_vkCreateDebugUtilsMessengerEXT) vkGetInstanceProcAddr(instance, "vkCreateDebugUtilsMessengerEXT");
        if (func == nullptr)
            return VK_ERROR_EXTENSION_NOT_PRESENT;
        else
            return func(instance, pCreateInfo, pAllocator, pDebugMessenger);
    }

    VkResult VkResourceManager::DestroyDebugUtilsMessengerEXT(
        VkInstance instance, 
        VkDebugUtilsMessengerEXT debugMessenger, 
        const VkAllocationCallbacks* pAllocator
    )
    {
        auto func = (PFN_vkDestroyDebugUtilsMessengerEXT)vkGetInstanceProcAddr(instance, "vkDestroyDebugUtilsMessengerEXT");
        if (func != nullptr)
        {
            func(instance, debugMessenger, pAllocator);
            return VK_SUCCESS;
        }
        else
            return VK_ERROR_EXTENSION_NOT_PRESENT;
    }

    void VkResourceManager::GetRequiredExtensions()
    {
        this->instanceExtensions.clear();

        UInt extensionCount = 0;
        vkEnumerateInstanceExtensionProperties(nullptr, &extensionCount, nullptr);

        std::vector<VkExtensionProperties> extensions(extensionCount);
        vkEnumerateInstanceExtensionProperties(nullptr, &extensionCount, extensions.data());

        UInt requiredExtensionCount = 0;
        const char** requiredExtensionNames = glfwGetRequiredInstanceExtensions(&requiredExtensionCount);

        for (UInt i = 0; i != requiredExtensionCount; ++i)
            this->instanceExtensions.push_back(requiredExtensionNames[i]);

        if (enabledValidation)
            instanceExtensions.push_back(VK_EXT_DEBUG_UTILS_EXTENSION_NAME);

        // 4189 for ATR_VERBOSE off
#pragma warning(disable: 4189)
        ATR_LOG_SECTION("Vulkan Available Extensions: (" + std::to_string(extensions.size()) + ")")
        for (const auto& extension : extensions)
            ATR_PRINT_VERBOSE(String("- ") + extension.extensionName)

        ATR_LOG_SECTION("Vulkan Required Extensions: (" + std::to_string(instanceExtensions.size()) + ")")
        for (const auto& extension : instanceExtensions)
            ATR_PRINT_VERBOSE(String("- ") + extension)
#pragma warning(default: 4189)

        // Check if all required extensions are available
        for (UInt i = 0; i != requiredExtensionCount; ++i)
            if (std::find_if(extensions.begin(), extensions.end(),
                [&](const VkExtensionProperties& ext) { return strcmp(ext.extensionName, requiredExtensionNames[i]) == 0; }
            ) == extensions.cend())
                throw Exception("Required extension not found: " + String(requiredExtensionNames[i]), ExceptionType::INIT_VULKAN);
    }

    void VkResourceManager::FindValidationLayers()
    {
        /// Validation Layers 
        if (enabledValidation)
        {
            UInt validatedLayerCount = 0;
            vkEnumerateInstanceLayerProperties(&validatedLayerCount, nullptr);

            std::vector<VkLayerProperties> availableValidatedLayers(validatedLayerCount);
            vkEnumerateInstanceLayerProperties(&validatedLayerCount, availableValidatedLayers.data());

            // 4189 for ATR_VERBOSE off
#pragma warning(disable: 4189)
            ATR_LOG_SECTION("Vulkan Available Layers: (" + std::to_string(validatedLayerCount) + ")")
            for (const auto& layer : availableValidatedLayers)
                ATR_PRINT_VERBOSE(String("- ") + layer.layerName)

            ATR_LOG_SECTION("Vulkan Required Layers: (" + std::to_string(this->validationLayers.size()) + ")")
            for (const auto& layer : this->validationLayers)
                ATR_PRINT_VERBOSE(String("- ") + layer)
#pragma warning(default: 4189)

            // Check if all validation layers are available
            for (auto& layerName : this->validationLayers)
                if (std::find_if(availableValidatedLayers.begin(), availableValidatedLayers.end(),
                    [&](const VkLayerProperties& layer) { return strcmp(layer.layerName, layerName) == 0; }
                ) == availableValidatedLayers.cend())
                    throw Exception("Validation layer " + String(layerName) + " not found", ExceptionType::INIT_VULKAN);
        }
    }

    Bool VkResourceManager::DeviceSuitable(VkPhysicalDevice device)
    {
        FindQueueFamilies(device);
        Bool deviceExtensionSupport = CheckDeviceExtensionSupport(device);
        QuerySwapChainSupport(device);

        Bool suitable =
            this->queueIndices.Complete() &&
            deviceExtensionSupport &&
            this->swapChainSupport.Adequate();

        if (suitable)
            return true;
        return false;
    }

    void VkResourceManager::FindQueueFamilies(VkPhysicalDevice device)
    {
        UInt queueFamilyCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, nullptr);

        std::vector<VkQueueFamilyProperties> queueFamilies(queueFamilyCount);
        vkGetPhysicalDeviceQueueFamilyProperties(device, &queueFamilyCount, queueFamilies.data());

        for (UInt i = 0; i != queueFamilies.size(); ++i)
        {
            if (queueFamilies[i].queueFlags & VK_QUEUE_GRAPHICS_BIT)
                this->queueIndices.indices[QueueFamilyIndices::GRAPHICS] = i;

            // if there exists a queue family that exclusively supports transfer operations, use it
            if (!(queueFamilies[i].queueFlags & VK_QUEUE_GRAPHICS_BIT) &&
                queueFamilies[i].queueFlags & VK_QUEUE_TRANSFER_BIT)
                this->queueIndices.indices[QueueFamilyIndices::TRANSFER] = i;

            VkBool32 presentSupport = false;
            vkGetPhysicalDeviceSurfaceSupportKHR(device, i, this->surface, &presentSupport);
            if (presentSupport)
                this->queueIndices.indices[QueueFamilyIndices::PRESENT] = i;

            if (this->queueIndices.Complete())
                break;
        }

        // Fallback: if there does not exist a queue exclusively for transfer operations, use the graphics queue
        if (this->queueIndices.indices[QueueFamilyIndices::TRANSFER] == std::nullopt)
            this->queueIndices.indices[QueueFamilyIndices::TRANSFER] = this->queueIndices.indices[QueueFamilyIndices::GRAPHICS];
    }

    // TODO replace repetitions of this code block
    void VkResourceManager::CreateImage(UInt width, UInt height, VkFormat format, VkImageTiling tiling, VkImageUsageFlags usage, VkMemoryPropertyFlags properties, VkImage& image, VkDeviceMemory& imageMemory)
    {
        VkImageCreateInfo imageInfo{};
        imageInfo.sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO;
        imageInfo.imageType = VK_IMAGE_TYPE_2D;
        imageInfo.extent.width = width;
        imageInfo.extent.height = height;
        imageInfo.extent.depth = 1;
        imageInfo.mipLevels = 1;
        imageInfo.arrayLayers = 1;
        imageInfo.format = format;
        imageInfo.tiling = tiling;
        imageInfo.initialLayout = VK_IMAGE_LAYOUT_UNDEFINED;
        imageInfo.usage = usage;
        imageInfo.samples = VK_SAMPLE_COUNT_1_BIT;
        imageInfo.sharingMode = VK_SHARING_MODE_EXCLUSIVE;

        if (vkCreateImage(device, &imageInfo, nullptr, &image) != VK_SUCCESS) {
            throw std::runtime_error("failed to create image!");
        }

        VkMemoryRequirements memRequirements;
        vkGetImageMemoryRequirements(device, image, &memRequirements);

        VkMemoryAllocateInfo allocInfo{};
        allocInfo.sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
        allocInfo.allocationSize = memRequirements.size;
        allocInfo.memoryTypeIndex = this->FindMemoryType(memRequirements.memoryTypeBits, properties);

        if (vkAllocateMemory(device, &allocInfo, nullptr, &imageMemory) != VK_SUCCESS) {
            throw std::runtime_error("failed to allocate image memory!");
        }

        vkBindImageMemory(device, image, imageMemory, 0);
    }

    VkImageView VkResourceManager::CreateImageView(VkImage image, VkFormat format, VkImageAspectFlags aspectFlags)
    {
        VkImageViewCreateInfo viewInfo{};
        viewInfo.sType = VK_STRUCTURE_TYPE_IMAGE_VIEW_CREATE_INFO;
        viewInfo.image = image;
        viewInfo.viewType = VK_IMAGE_VIEW_TYPE_2D;
        viewInfo.format = format;
        viewInfo.subresourceRange.aspectMask = aspectFlags;
        viewInfo.subresourceRange.baseMipLevel = 0;
        viewInfo.subresourceRange.levelCount = 1;
        viewInfo.subresourceRange.baseArrayLayer = 0;
        viewInfo.subresourceRange.layerCount = 1;

        VkImageView imageView;
        if (vkCreateImageView(device, &viewInfo, nullptr, &imageView) != VK_SUCCESS) {
            throw std::runtime_error("failed to create texture image view!");
        }

        return imageView;
    }

    VkFormat VkResourceManager::FindSupportedFormat(const std::vector<VkFormat>& candidates, VkImageTiling tiling, VkFormatFeatureFlags features)
    {
        for (VkFormat format : candidates)
        {
            VkFormatProperties props;
            vkGetPhysicalDeviceFormatProperties(this->physicalDevice, format, &props);

            if (tiling == VK_IMAGE_TILING_LINEAR && (features & props.linearTilingFeatures) == features)
                return format;
            else if (tiling == VK_IMAGE_TILING_OPTIMAL && (features & props.optimalTilingFeatures) == features)
                return format;
        }

        throw Exception("Failed to find supported format", ExceptionType::INIT_PIPELINE);
    }

    Bool VkResourceManager::CheckDeviceExtensionSupport(VkPhysicalDevice device)
    {
        UInt deviceExtensionCount = 0;
        vkEnumerateDeviceExtensionProperties(device, nullptr, &deviceExtensionCount, nullptr);

        std::vector<VkExtensionProperties> availableDeviceExtensions(deviceExtensionCount);
        vkEnumerateDeviceExtensionProperties(device, nullptr, &deviceExtensionCount, availableDeviceExtensions.data());

        for (auto& extensionName : this->deviceExtensions)
        {
            if (std::find_if(availableDeviceExtensions.begin(), availableDeviceExtensions.end(),
                [&](const VkExtensionProperties& ext) { return strcmp(ext.extensionName, extensionName) == 0; }
            ) == availableDeviceExtensions.cend())
                return false;
        }

        return true;
    }

    void VkResourceManager::QuerySwapChainSupport(VkPhysicalDevice device)
    {
        // Capabilities
        vkGetPhysicalDeviceSurfaceCapabilitiesKHR(device, this->surface, &this->swapChainSupport.capabilities);

        // Formats
        UInt formatCount = 0;
        vkGetPhysicalDeviceSurfaceFormatsKHR(device, this->surface, &formatCount, nullptr);

        if (formatCount != 0)
        {
            this->swapChainSupport.formats.resize(formatCount);
            vkGetPhysicalDeviceSurfaceFormatsKHR(device, this->surface, &formatCount, this->swapChainSupport.formats.data());
        }

        // Present Modes
        UInt presentModeCount = 0;
        vkGetPhysicalDeviceSurfacePresentModesKHR(device, this->surface, &presentModeCount, nullptr);

        if (presentModeCount != 0)
        {
            this->swapChainSupport.presentModes.resize(presentModeCount);
            vkGetPhysicalDeviceSurfacePresentModesKHR(device, this->surface, &presentModeCount, this->swapChainSupport.presentModes.data());
        }
    }

    void VkResourceManager::ConfigureSwapChain(SwapChainSupportDetails support)
    {
        // Choose Swap Surface Format
        Bool found = false;
        for (const auto& format : support.formats)
        {
            if (format.format == VkResourceManager::defaultColorFormat && format.colorSpace == VK_COLOR_SPACE_SRGB_NONLINEAR_KHR)
            {
                this->swapChainConfig.format = format;
                found = true;
                break;
            }
        }

        if (!found)
            this->swapChainConfig.format = support.formats[0];

        // Choose Swapchain Presentation Mode
        found = false;
        for (const auto& mode : support.presentModes)
        {
            if (mode == VK_PRESENT_MODE_MAILBOX_KHR)
            {
                this->swapChainConfig.presentMode = mode;
                found = true;
                break;
            }
        }

        if (!found)
            this->swapChainConfig.presentMode = VK_PRESENT_MODE_FIFO_KHR;

        // Choose Swap Extent
        if (support.capabilities.currentExtent.width != std::numeric_limits<UInt>::max())
            this->swapChainConfig.extent = support.capabilities.currentExtent;
        else
        {
            Int width, height;
            glfwGetFramebufferSize(this->window, &width, &height);

            VkExtent2D actualExtent = {
                .width = static_cast<UInt>(width),
                .height = static_cast<UInt>(height)
            };

            actualExtent.width = std::max(support.capabilities.minImageExtent.width,
                std::min(support.capabilities.maxImageExtent.width, actualExtent.width));
            actualExtent.height = std::max(support.capabilities.minImageExtent.height,
                std::min(support.capabilities.maxImageExtent.height, actualExtent.height));

            this->swapChainConfig.extent = actualExtent;
        }
    }

    void VkResourceManager::CreateBuffer(VkDeviceSize size, VkBufferUsageFlags usage, VkMemoryPropertyFlags properties, VkBuffer& buffer, VkDeviceMemory& bufferMemory)
    {
        VkBufferCreateInfo bufferInfo = {
            .sType = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO,
            .size = size,
            .usage = usage,
            .sharingMode = VK_SHARING_MODE_EXCLUSIVE
        };

        if (vkCreateBuffer(this->device, &bufferInfo, nullptr, &buffer) != VK_SUCCESS)
            throw Exception("Failed to create buffer", ExceptionType::INIT_PIPELINE);

        VkMemoryRequirements memRequirements;
        vkGetBufferMemoryRequirements(this->device, buffer, &memRequirements);

        VkMemoryAllocateInfo allocInfo = {
            .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
            .allocationSize = memRequirements.size,
            .memoryTypeIndex = this->FindMemoryType(memRequirements.memoryTypeBits, properties)
        };

        // TODO implement a custom allocator for handling multiple memory allocation, with different offsets
        //   refer to https://vulkan-tutorial.com/Vertex_buffers/Staging_buffer#page_Conclusion
        //   Also create index/vertex buffers in a single allocation
        //   refer to https://developer.nvidia.com/vulkan-memory-management
        VkResult result = vkAllocateMemory(this->device, &allocInfo, nullptr, &bufferMemory);
        if (result != VK_SUCCESS)
            throw Exception("Failed to allocate buffer memory", ExceptionType::INIT_BUFFER);

        vkBindBufferMemory(this->device, buffer, bufferMemory, 0);
    }

    void VkResourceManager::CreateStagingBuffer(VkDeviceSize size, VkBuffer& stagingBuffer, VkDeviceMemory& stagingBufferMemory)
    {
        this->CreateBuffer(
            size,
            VK_BUFFER_USAGE_TRANSFER_SRC_BIT,           // Will be transferred to the actual vertex buffer on device
            VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
            // Host coherency implies that changes to the memory on the host will NOTIFY (but not necessarily immediately visible to) the device
            // Device visibility depends on its actual memory model and is ensured by the driver
            stagingBuffer,
            stagingBufferMemory
        );
    }

    void VkResourceManager::CopyBuffer(VkBuffer src, VkBuffer dst, VkDeviceSize size)
    {
        VkCommandBufferAllocateInfo allocInfo = {
            .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
            .commandPool = this->GetTransferCommandPool(),
            .level = VK_COMMAND_BUFFER_LEVEL_PRIMARY,
            .commandBufferCount = 1
        };

        if (!this->queueIndices.SeparateTransferQueue())
            allocInfo.commandPool = this->graphicsCommandPool;

        VkCommandBuffer copyCommandBuffer;
        vkAllocateCommandBuffers(this->device, &allocInfo, &copyCommandBuffer);

        VkCommandBufferBeginInfo beginInfo = {
            .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
            .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT
        };

        // Recording command buffer for transfer queue
        if (vkBeginCommandBuffer(copyCommandBuffer, &beginInfo) != VK_SUCCESS)
            throw Exception("Failed to begin recording command buffer for transfer", ExceptionType::UPDATE_MEMORY);
            
            VkBufferCopy copyRegion = {
                .srcOffset = 0,
                .dstOffset = 0,
                .size = size
            };
            vkCmdCopyBuffer(copyCommandBuffer, src, dst, 1, &copyRegion);

        if (vkEndCommandBuffer(copyCommandBuffer) != VK_SUCCESS)
            throw Exception("Failed to end recording command buffer for transfer", ExceptionType::UPDATE_MEMORY);

        VkSubmitInfo submitInfo = {
            .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
            .commandBufferCount = 1,
            .pCommandBuffers = &copyCommandBuffer
        };

        vkQueueSubmit(this->queues[QueueFamilyIndices::TRANSFER], 1, &submitInfo, VK_NULL_HANDLE);
        vkQueueWaitIdle(this->queues[QueueFamilyIndices::TRANSFER]);

        vkFreeCommandBuffers(this->device, this->GetTransferCommandPool(), 1, &copyCommandBuffer);
    }

    void VkResourceManager::RetrieveSwapChainImages()
    {
        UInt imageCount = 0;
        vkGetSwapchainImagesKHR(this->device, this->swapchain, &imageCount, nullptr);

        this->swapchainImages.resize(imageCount);
        vkGetSwapchainImagesKHR(this->device, this->swapchain, &imageCount, this->swapchainImages.data());

        ATR_PRINT_VERBOSE("Retrieved " + std::to_string(imageCount) + " images in total in the swapchain.")
    }

    inline VkFormat VkResourceManager::FindDepthFormat()
    {
        return FindSupportedFormat(
            { VK_FORMAT_D32_SFLOAT, VK_FORMAT_D32_SFLOAT_S8_UINT, VK_FORMAT_D24_UNORM_S8_UINT },
            VK_IMAGE_TILING_OPTIMAL,
            VK_FORMAT_FEATURE_DEPTH_STENCIL_ATTACHMENT_BIT
        );
    }

    inline VkCommandPool& VkResourceManager::GetTransferCommandPool()
    {
        if (this->queueIndices.SeparateTransferQueue())
            return this->transferCommandPool;
        else
            return this->graphicsCommandPool;
    }

    inline bool VkResourceManager::HasStencilComponent(VkFormat format)
    {
        return format == VK_FORMAT_D32_SFLOAT_S8_UINT || format == VK_FORMAT_D24_UNORM_S8_UINT;
    }

    VkCommandBuffer VkResourceManager::CreateCommandBuffer(VkCommandBufferLevel level, VkCommandPool commandPool, Bool startRecording)
    {
        VkCommandBufferAllocateInfo allocInfo = 
        {
            .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO,
            .commandPool = commandPool,
            .level = level,
            .commandBufferCount = 1
        };
        VkCommandBuffer commandBuffer;

        if (vkAllocateCommandBuffers(this->device, &allocInfo, &commandBuffer) != VK_SUCCESS)
            throw Exception("Failed to allocate command buffer", ExceptionType::INIT_BUFFER);

        if (startRecording)
        {
            VkCommandBufferBeginInfo beginInfo =
            {
                .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
                .flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT
            };

            if (vkBeginCommandBuffer(commandBuffer, &beginInfo) != VK_SUCCESS)
                throw Exception("Failed to begin recording command buffer", ExceptionType::INIT_BUFFER);
        }

        return commandBuffer;
    }

    void VkResourceManager::FlushCommandBuffer(VkCommandBuffer commandBuffer, VkCommandPool commandPool, VkQueue queue, Bool free)
    {
        if (commandBuffer == VK_NULL_HANDLE)
            return; 

        if (vkEndCommandBuffer(commandBuffer) != VK_SUCCESS)
            throw Exception("Failed to end recording command buffer", ExceptionType::UPDATE_MEMORY);

        VkSubmitInfo submitInfo = 
        {
            .sType = VK_STRUCTURE_TYPE_SUBMIT_INFO,
            .commandBufferCount = 1,
            .pCommandBuffers = &commandBuffer
        };

        VkFenceCreateInfo fenceInfo = 
        {
            .sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO,
            .flags = 0
        };
        VkFence fence;

        if (vkCreateFence(this->device, &fenceInfo, nullptr, &fence) != VK_SUCCESS)
            throw Exception("Failed to create fence", ExceptionType::UPDATE_MEMORY);

        if (vkQueueSubmit(queue, 1, &submitInfo, fence) != VK_SUCCESS)
            throw Exception("Failed to submit command buffer", ExceptionType::UPDATE_MEMORY);

        if (vkWaitForFences(this->device, 1, &fence, VK_TRUE, UINT64_MAX) != VK_SUCCESS)
            throw Exception("Failed to wait for fence", ExceptionType::UPDATE_MEMORY);

        vkDestroyFence(this->device, fence, nullptr);

        if (free)
            vkFreeCommandBuffers(this->device, commandPool, 1, &commandBuffer);
    }

    void VkResourceManager::InsertImageMemoryBarrier(
        VkCommandBuffer commandBuffer, 
        VkImage image, 
        VkAccessFlags srcAccessMask, 
        VkAccessFlags dstAccessMask, 
        VkImageLayout oldLayout, 
        VkImageLayout newLayout, 
        VkPipelineStageFlags srcStageMask, 
        VkPipelineStageFlags dstStageMask, 
        VkImageSubresourceRange subresourceRange)
    {
        VkImageMemoryBarrier barrier = 
        {
            .sType = VK_STRUCTURE_TYPE_IMAGE_MEMORY_BARRIER,
            .srcAccessMask = srcAccessMask,
            .dstAccessMask = dstAccessMask,
            .oldLayout = oldLayout,
            .newLayout = newLayout,
            .srcQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED,
            .dstQueueFamilyIndex = VK_QUEUE_FAMILY_IGNORED,
            .image = image,
            .subresourceRange = subresourceRange
        };

        vkCmdPipelineBarrier
        (
            commandBuffer,
            srcStageMask,
            dstStageMask,
            0,
            0, nullptr,
            0, nullptr,
            1, &barrier
        );
    }

    void VkResourceManager::RecordCommandBuffer(VkCommandBuffer commandBuffer, UInt imageIndex)
    {
        VkCommandBufferBeginInfo beginInfo = {
            .sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO,
            .flags = 0,
            .pInheritanceInfo = nullptr
        };

        // Recording command for rendering
        if (vkBeginCommandBuffer(commandBuffer, &beginInfo) != VK_SUCCESS)
            throw Exception("Failed to begin recording command buffer for rendering", ExceptionType::INIT_PIPELINE);

        std::array<VkClearValue, 2> clearValues = { VkResourceManager::defaultClearValue, VkResourceManager::defaultDepthClearValue };

        VkRenderPassBeginInfo renderPassInfo = {
            .sType = VK_STRUCTURE_TYPE_RENDER_PASS_BEGIN_INFO,
            .renderPass = this->renderPass,
            .framebuffer = this->swapchainFrameBuffers[imageIndex],
            .renderArea = {
                .offset = { 0, 0 },
                .extent = this->swapChainConfig.extent
            },
            .clearValueCount = static_cast<UInt>(clearValues.size()),
            .pClearValues = clearValues.data()
        };

        vkCmdBeginRenderPass(commandBuffer, &renderPassInfo, VK_SUBPASS_CONTENTS_INLINE);

            vkCmdBindPipeline(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, this->graphicsPipeline);

            VkViewport viewport;
            viewport.x = 0;
            viewport.y = 0;
            viewport.width = (Float)this->swapChainConfig.extent.width;
            viewport.height = (Float)this->swapChainConfig.extent.height;
            viewport.minDepth = 0.0f;
            viewport.maxDepth = 1.0f;
            vkCmdSetViewport(this->graphicsCommandBuffers[this->currentFrameIndex], 0, 1, &viewport);

            VkRect2D scissor;
            scissor.offset = { 0, 0 };
            scissor.extent = this->swapChainConfig.extent;
            vkCmdSetScissor(this->graphicsCommandBuffers[this->currentFrameIndex], 0, 1, &scissor);

            VkBuffer vertexBuffers[] = { this->vertexBuffer };
            VkDeviceSize offsets[] = { 0 };
            vkCmdBindVertexBuffers(commandBuffer, 0, 1, vertexBuffers, offsets);
            vkCmdBindIndexBuffer(commandBuffer, this->indexBuffer, 0, VK_INDEX_TYPE_UINT32);

            vkCmdBindDescriptorSets(commandBuffer, VK_PIPELINE_BIND_POINT_GRAPHICS, this->pipelineLayout, 0, 1, &this->descriptorSets[this->currentFrameIndex], 0, nullptr);

            vkCmdDrawIndexed(commandBuffer, static_cast<UInt>(this->mesh.GetIndices().size()), 1, 0, 0, 0);

        vkCmdEndRenderPass(commandBuffer);

        if (vkEndCommandBuffer(commandBuffer) != VK_SUCCESS)
            throw Exception("Failed to end recording command buffer for rendering", ExceptionType::INIT_PIPELINE);
    }

    UInt VkResourceManager::FindMemoryType(UInt typeFilter, VkMemoryPropertyFlags properties)
    {
        VkPhysicalDeviceMemoryProperties memProperties;
        vkGetPhysicalDeviceMemoryProperties(this->physicalDevice, &memProperties);

        for (UInt i = 0; i != memProperties.memoryTypeCount; ++i)
            if (
                typeFilter & (1 << i) &&                                                    // Suitable for the buffer
                (memProperties.memoryTypes[i].propertyFlags & properties) == properties     // Suitable for data
            )
                return i;

        throw Exception("Failed to find suitable memory type", ExceptionType::UPDATE_MEMORY);
    }

    void VkResourceManager::UpdateUniformBuffer(UInt currentFrameIndex)
    {
        static auto startTime = std::chrono::high_resolution_clock::now();

        UniformBufferObject ubo;
        //ubo.model = glm::rotate(glm::mat4(1.0f), time * glm::radians(90.f), glm::vec3(0.0f, 0.0f, 1.0f));
        ubo.model = glm::mat4(1.0f);
        ubo.view = glm::lookAt(glm::vec3(2.f, 2.f, 2.f), glm::vec3(0.f, 0.f, 0.f), glm::vec3(0.f, 0.f, 1.f));
        auto swapExtent = this->swapChainConfig.extent;
        ubo.proj = glm::perspective(glm::radians(45.f), static_cast<float>(swapExtent.width) / static_cast<float>(swapExtent.height), 0.1f, 10.f);

        ubo.proj[1][1] *= -1;       // Vulkan designates the origin of an image to be the upper-left vertex
        
        // TODO learn about "push constants" for improving efficiency
        memcpy(this->uniformBufferMappedMemory[currentFrameIndex], &ubo, sizeof(ubo));
    } 

    void VkResourceManager::RecreateSwapchain()
    {
        // If window is minimized, pause update
        int width = 0, height = 0;
        glfwGetFramebufferSize(this->window, &width, &height);
        while (width == 0 || height == 0)
        {
            // TODO get rid of busy waiting here
            glfwGetFramebufferSize(this->window, &width, &height);
            glfwWaitEvents();
        }

        vkDeviceWaitIdle(this->device);

        CleanUpSwapchain();

        CreateSwapchain();
        CreateImageViews();
        CreateDepthBuffer();
        CreateFrameBuffers();
    }

    void VkResourceManager::UpdateImageBuffers()
    {
        ATR_LOG("Updating Mesh Infos...")
        VkDeviceSize indexBufferSize = sizeof(UInt) * this->mesh.GetIndices().size();

        VkBuffer indexStagingBuffer;
        VkDeviceMemory indexStagingBufferMemory;

        this->CreateStagingBuffer(indexBufferSize, indexStagingBuffer, indexStagingBufferMemory);

        void* indexData;
        vkMapMemory(this->device, indexStagingBufferMemory, 0, indexBufferSize, 0, &indexData);
        memcpy(indexData, this->mesh.GetIndices().data(), static_cast<size_t>(indexBufferSize));
        vkUnmapMemory(this->device, indexStagingBufferMemory);
        
        this->CopyBuffer(indexStagingBuffer, this->indexBuffer, indexBufferSize);

        vkDestroyBuffer(this->device, indexStagingBuffer, nullptr);
        vkFreeMemory(this->device, indexStagingBufferMemory, nullptr);

        VkDeviceSize vertexBufferSize = sizeof(Vertex) * this->mesh.GetVertices().size();

        VkBuffer vertexStagingBuffer;
        VkDeviceMemory vertexStagingBufferMemory;

        this->CreateStagingBuffer(vertexBufferSize, vertexStagingBuffer, vertexStagingBufferMemory);

        void* vertexData;
        vkMapMemory(this->device, vertexStagingBufferMemory, 0, vertexBufferSize, 0, &vertexData);
        memcpy(vertexData, this->mesh.GetVertices().data(), static_cast<size_t>(vertexBufferSize));
        vkUnmapMemory(this->device, vertexStagingBufferMemory);

        this->CopyBuffer(vertexStagingBuffer, this->vertexBuffer, vertexBufferSize);

        vkDestroyBuffer(this->device, vertexStagingBuffer, nullptr);
        vkFreeMemory(this->device, vertexStagingBufferMemory, nullptr);
    }

    void VkResourceManager::CompileShaders()
    {
        const String& path = this->relLocation;
        const String sep = (path == "") ? "" : "\\";
        const String compilerPath = path + sep + "..\\vendor\\bin\\Windows\\glslc.exe";

        ATR::OS::Execute("rmdir /s /q bin\\shaders");
        ATR::OS::Execute("mkdir bin");
        ATR::OS::Execute("mkdir bin\\shaders");
        ATR::OS::Execute(compilerPath + " " + "..\\B-Spline-IPC\\shaders\\shader.vert -o bin\\shaders\\vert.spv");
        ATR::OS::Execute(compilerPath + " " + "..\\B-Spline-IPC\\shaders\\shader.frag -o bin\\shaders\\frag.spv");
    }

    std::vector<char> VkResourceManager::ReadShaderCode(const char* path)
    {
        std::ifstream file(path, std::ios::ate | std::ios::binary);
        if (!file.is_open())
            throw Exception("Failed to open file: " + String(path), ExceptionType::INIT_SHADER);

        size_t fileSize = (size_t)file.tellg();
        std::vector<char> buffer(fileSize);
        file.seekg(0);
        file.read(buffer.data(), fileSize);
        file.close();

        return buffer;
    }

    VkShaderModule VkResourceManager::CreateShaderModule(const std::vector<char>& code)
    {
        VkShaderModuleCreateInfo createInfo{
            .sType = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO,
            .codeSize = code.size(),
            .pCode = reinterpret_cast<const UInt*>(code.data())
        };

        VkShaderModule shaderModule;
        if (vkCreateShaderModule(this->device, &createInfo, nullptr, &shaderModule) != VK_SUCCESS)
            throw Exception("Failed to create shader module", ExceptionType::INIT_SHADER);

        return shaderModule;
    }

    void VkResourceManager::SaveScreenshot(String filename, UInt imageIndex)
    {
        bool supportsBlit = true;

        // Check blit support for source and destination
        VkFormatProperties formatProps;

        // Check if the device supports blitting from optimal images (the swapchain images are in optimal format)
        vkGetPhysicalDeviceFormatProperties(physicalDevice, VkResourceManager::defaultColorFormat, &formatProps);
        if (!(formatProps.optimalTilingFeatures & VK_FORMAT_FEATURE_BLIT_SRC_BIT)) {
            std::cerr << "Device does not support blitting from optimal tiled images, using copy instead of blit!" << std::endl;
            supportsBlit = false;
        }

        // Check if the device supports blitting to linear images
        vkGetPhysicalDeviceFormatProperties(physicalDevice, VK_FORMAT_R8G8B8A8_SRGB, &formatProps);
        if (!(formatProps.linearTilingFeatures & VK_FORMAT_FEATURE_BLIT_DST_BIT)) {
            std::cerr << "Device does not support blitting to linear tiled images, using copy instead of blit!" << std::endl;
            supportsBlit = false;
        }

        // Duplicate source Image
        VkImage srcImage = this->swapchainImages[imageIndex];

        // Prepare image pending copy
        VkImageCreateInfo imageCreateInfo = 
        {
            .sType = VK_STRUCTURE_TYPE_IMAGE_CREATE_INFO,
            .imageType = VK_IMAGE_TYPE_2D,
            .format = VK_FORMAT_R8G8B8A8_SRGB,
            .extent = 
            { 
                .width = this->swapChainConfig.extent.width, 
                .height = this->swapChainConfig.extent.height, 
                .depth = 1 
            },
            .mipLevels = 1,
            .arrayLayers = 1,
            .samples = VK_SAMPLE_COUNT_1_BIT,
            .tiling = VK_IMAGE_TILING_LINEAR,
            .usage = VK_IMAGE_USAGE_TRANSFER_DST_BIT,
            .sharingMode = VK_SHARING_MODE_EXCLUSIVE,
            .initialLayout = VK_IMAGE_LAYOUT_UNDEFINED
        };

        VkImage dstImage;
        if (vkCreateImage(this->device, &imageCreateInfo, nullptr, &dstImage) != VK_SUCCESS)
            throw Exception("Failed to create image for screenshot", ExceptionType::UPDATE_MEMORY);

        // Create memory for the image
        VkMemoryRequirements memRequirements;
        vkGetImageMemoryRequirements(this->device, dstImage, &memRequirements);

        VkMemoryAllocateInfo allocInfo = 
        {
            .sType = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO,
            .allocationSize = memRequirements.size,
            .memoryTypeIndex = this->FindMemoryType(memRequirements.memoryTypeBits, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT)
        };

        VkDeviceMemory dstImageMemory;

        if (vkAllocateMemory(this->device, &allocInfo, nullptr, &dstImageMemory) != VK_SUCCESS)
            throw Exception("Failed to allocate memory for screenshot image", ExceptionType::UPDATE_MEMORY);

        if (vkBindImageMemory(this->device, dstImage, dstImageMemory, 0) != VK_SUCCESS)
            throw Exception("Failed to bind memory to screenshot image", ExceptionType::UPDATE_MEMORY);

        // Prepare command buffer for image copy
        VkCommandBuffer copyCmdBuffer = this->CreateCommandBuffer(VK_COMMAND_BUFFER_LEVEL_PRIMARY, this->graphicsCommandPool, true);

        /// Transition the images to correct layout for transfer
        // Transition source image to transfer source layout
        this->InsertImageMemoryBarrier(
            copyCmdBuffer,
            srcImage,
            VK_ACCESS_MEMORY_READ_BIT,
            VK_ACCESS_TRANSFER_READ_BIT,
            VK_IMAGE_LAYOUT_PRESENT_SRC_KHR, 
            VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
            VK_PIPELINE_STAGE_TRANSFER_BIT, 
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VkImageSubresourceRange{ VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 }
        );

        // Transition destination image to transfer destination layout
        this->InsertImageMemoryBarrier(
            copyCmdBuffer,
            dstImage,
            0,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_IMAGE_LAYOUT_UNDEFINED,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VkImageSubresourceRange{ VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 }
        );

        // Conduct the actual transfer based on whether blit is supported
        if (supportsBlit)
        {
            VkOffset3D blitSize =
            {
                .x = static_cast<int32_t>(this->swapChainConfig.extent.width),
                .y = static_cast<int32_t>(this->swapChainConfig.extent.height),
                .z = 1
            };

            VkImageBlit imageBlitRegion =
            {
                .srcSubresource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 },
                .srcOffsets = { { 0, 0, 0 }, blitSize },
                .dstSubresource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 },
                .dstOffsets = { { 0, 0, 0 }, blitSize }
            };

            vkCmdBlitImage(
                copyCmdBuffer,
                srcImage, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                dstImage, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                1, &imageBlitRegion,
                VK_FILTER_NEAREST
            );
        }
        else
        {
            VkImageCopy imageCopyRegion = 
            {
                .srcSubresource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 },
                .srcOffset = { 0, 0, 0 },
                .dstSubresource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0, 1 },
                .dstOffset = { 0, 0, 0 },
                .extent = { this->swapChainConfig.extent.width, this->swapChainConfig.extent.height, 1 }
            };

            vkCmdCopyImage(
                copyCmdBuffer,
                srcImage, VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
                dstImage, VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
                1, &imageCopyRegion
            );
        }

        /// After doing the actual copy, transfer back the images to their original layout
        // Transition back the source image to present layout
        this->InsertImageMemoryBarrier(
            copyCmdBuffer,
            srcImage,
            VK_ACCESS_TRANSFER_READ_BIT,
            VK_ACCESS_MEMORY_READ_BIT,
            VK_IMAGE_LAYOUT_TRANSFER_SRC_OPTIMAL,
            VK_IMAGE_LAYOUT_PRESENT_SRC_KHR,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VkImageSubresourceRange{ VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 }
        );

        // Transition the destination image to general layout for host access
        this->InsertImageMemoryBarrier(
            copyCmdBuffer,
            dstImage,
            VK_ACCESS_TRANSFER_WRITE_BIT,
            VK_ACCESS_MEMORY_READ_BIT,
            VK_IMAGE_LAYOUT_TRANSFER_DST_OPTIMAL,
            VK_IMAGE_LAYOUT_GENERAL,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VK_PIPELINE_STAGE_TRANSFER_BIT,
            VkImageSubresourceRange{ VK_IMAGE_ASPECT_COLOR_BIT, 0, 1, 0, 1 }
        );

        // Flush command buffer
        this->FlushCommandBuffer(copyCmdBuffer, this->graphicsCommandPool, this->queues[QueueFamilyIndices::GRAPHICS], true);

        // Get Layout of the image
        VkImageSubresource subResource = { VK_IMAGE_ASPECT_COLOR_BIT, 0, 0 };
        VkSubresourceLayout subResourceLayout;
        vkGetImageSubresourceLayout(this->device, dstImage, &subResource, &subResourceLayout);

        // Map the image memory to CPU accessible memory
        const char* data;
        vkMapMemory(this->device, dstImageMemory, 0, VK_WHOLE_SIZE, 0, (void**)&data);
        data += subResourceLayout.offset;

        // Write to png
        Int info;
        info = stbi_write_png(
            filename.c_str(),
            this->swapChainConfig.extent.width, this->swapChainConfig.extent.height, 
            4,
            data,
            static_cast<Int>(subResourceLayout.rowPitch)
        );
    }

    // Vulkan callback, can allow unused parameters
#pragma warning(disable: 4100)
    void VkResourceManager::FramebufferResizeCallback(GLFWwindow* window, int width, int height)
    {
        auto app = reinterpret_cast<VkResourceManager*>(glfwGetWindowUserPointer(window));
        ATR_LOG_VERBOSE("Window Resizing... new width: " << width << ", height: " << height)
        app->frameBufferResized = true;
    }
#pragma warning(default: 4100)
}