#pragma once
#include "atrfwd.h"

#include <vector>

namespace ATR
{
    struct SwapChainSupportDetails
    {
        // Contains information about compatibility of swap chain and the surface

        VkSurfaceCapabilitiesKHR capabilities;                  // Surface capabilities: images in swap chain; min/max width/height
        std::vector<VkSurfaceFormatKHR> formats;                // Surface formats: pixel format, color space
        std::vector<VkPresentModeKHR> presentModes;             // Presentation modes: how image swapping is implemented

        inline Bool Adequate() { return !formats.empty() && !presentModes.empty(); }

        // Too verbose for printing... not defining ostream operator
    };
}
