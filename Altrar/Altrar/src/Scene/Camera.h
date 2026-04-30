#pragma once

#include "atrfwd.h"

namespace ATR
{
    struct Camera
    {
        Vec3 pos;
        Vec3 lookAt;
        Vec3 up;

        Float width;
        Float height;

        Float nearClip;
        Float farClip;
    };

}