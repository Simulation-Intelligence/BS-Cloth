#pragma once
#include "atrfwd.h"

namespace ATR
{

#define ATR_UNIFORM_VEC2 alignas(8)  Vec2
#define ATR_UNIFORM_VEC3 alignas(16) Vec3
#define ATR_UNIFORM_VEC4 alignas(16) Vec4
#define ATR_UNIFORM_MAT4 alignas(16) Mat4

    struct UniformBufferObject
    {
        ATR_UNIFORM_MAT4 model;
        ATR_UNIFORM_MAT4 view;
        ATR_UNIFORM_MAT4 proj;
    };
}
