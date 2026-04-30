#pragma once

#include "bsfwd.h"

namespace BSIPC
{
    inline Bool InClosedInterval(Float x, Float a, Float b)
    {
        BSIPC_ASSERT(a <= b, "Invalid interval")
        return x >= a && x <= b;
    }

    inline Bool InOpenInterval(Float x, Float a, Float b)
    {
        BSIPC_ASSERT(a < b, "Invalid interval")
        return x > a && x < b;
    }

    inline Float trigArea(Vec2 v1, Vec2 v2, Vec2 v3)
    {
        Float area = 0.5 * std::abs(
            v1[0] * (v2[1] - v3[1]) +
            v2[0] * (v3[1] - v1[1]) +
            v3[0] * (v1[1] - v2[1])
        );
        return area;
    }

    inline LUInt PairHashKey(UInt x, UInt y)
    {
        return (static_cast<LUInt>(x) << 32) | y;
    }

    inline LInt PairHashKey(Int x, Int y)
    {
        return (static_cast<LInt>(x) << 32) | y;
    }
}