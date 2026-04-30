#pragma once

#include "stdint.h"
#include <string>

#define GLM_FORCE_RADIANS
#define GLM_FORCE_DEPTH_ZERO_TO_ONE
#define GLM_FORCE_DEFAULT_ALIGNED_GENTYPES
#include "glm/glm.hpp"

namespace ATR
{
    // STL Types
    using Int = int32_t;
    using UInt = uint32_t;
    using LInt = int64_t;
    using LUInt = uint64_t;
    using Char = char8_t;
    using Bool = bool;
    using String = std::string;
    using Float = float;

    // Common GLM Types
    using Vec2 = glm::vec2;
    using Vec3 = glm::vec3;
    using Vec4 = glm::vec4;
    using Mat2 = glm::mat2x2;
    using Mat3 = glm::mat3x3;
    using Mat4 = glm::mat4x4;
}