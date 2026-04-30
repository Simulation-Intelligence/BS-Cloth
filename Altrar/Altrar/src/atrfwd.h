#pragma once

/// Forward declaration of constructions needed by headers, and should enter into the precompiled header
// Can only include files in src/Core/
// Headers in src/Core should not include atrfwd.h, and should instead include the corresponding header file directly

// Types
// GLM is included here
#include "ATRType.h"

// OS Specifics and Outputs
#include "ATRInfo.h"

// Global Defines
#if defined ATR_DEBUG


#elif defined ATR_RELEASE


#endif 

// STL Containers
#include <vector>
#include <array>
#include <set>

// STL Descriptives
#include <optional>

// STL Tools
#include <algorithm>

// Vulkan 
#include "vulkan/vulkan.h"
