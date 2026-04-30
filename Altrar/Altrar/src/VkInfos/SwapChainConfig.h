#pragma once
#include "atrfwd.h"

namespace ATR
{
    struct SwapChainConfig
    {
        VkSurfaceFormatKHR format;
        VkPresentModeKHR presentMode;
        VkExtent2D extent;
    };
}
