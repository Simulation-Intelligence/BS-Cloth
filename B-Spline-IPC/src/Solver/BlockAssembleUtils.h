#pragma once

#include "bsfwd.h"

namespace BSIPC
{
    enum class LinearHessBlockType
    {
        Block3x3 = 0,
        Block6x6 = 1,
        Block9x9 = 2,
        Block12x12 = 3
    };

    struct LinearHessBlockSource
    {
        LinearHessBlockType blockType;
        UInt vecIndex;
        IVec2 blockIndices;
    };

    struct BSHessBlockSource
    {
        const Mat<3, 3>* sourceBlock;
        Float coeff;
    };

    struct ElasticityGradBlockSource
    {
        UInt globalPatchIdx;
        UInt blockRowIdx;
    };

    struct ElasticityHessBlockSource
    {
        UInt globalPatchIdx;
        UInt blockRowIdx;
        UInt blockColIdx;
    };

    struct ElasticityHessEntrySource
    {
        UInt globalPatchIdx;
        UInt rowIdx;
        UInt colIdx;
    };
}
