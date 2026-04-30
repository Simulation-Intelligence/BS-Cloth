#pragma once
#include "bsfwd.h"

#include "Geometry/BSpline/BSSurface.h"

namespace BSIPC
{
    struct BSTargetInfo
    {
        BSTargetInfo() = default;
        BSTargetInfo(BSSurface* target, UInt index) :
            target(target), index(index),
            contactTrigVertIndexOffset(std::nullopt), trigVertCnt(std::nullopt),
            renderTrigVertIndexOffset(std::nullopt),
            bsVertIndexOffset(std::nullopt), bsPatchIndexOffset(std::nullopt)
        {  }

        BSSurface* const target;
        const UInt index;
        std::optional<UInt> contactTrigVertIndexOffset;             // The start index of target in Solver::contactMesh
        std::optional<UInt> trigVertCnt;                            // The number of vertices in target in Solver::contactMesh
        std::optional<UInt> renderTrigVertIndexOffset;              // The start index of target in Solver::renderMesh
        std::optional<UInt> bsVertIndexOffset;                      // The start index of target in sequenced BS surfaces
        std::optional<UInt> bsPatchIndexOffset;                     // The start index of patch in sequenced BS surfaces
        // no need to have bsVertCnt, as it can be calculated from BSSurface::GetControlPointCnt()

        inline Bool ValidateTarget() const { return target != nullptr; }
        inline Bool ValidateTrig() const { return contactTrigVertIndexOffset.has_value() && trigVertCnt.has_value(); }
        inline Bool Validate() const { return ValidateTarget() && ValidateTrig(); }
    };
}
