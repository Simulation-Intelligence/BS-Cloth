#pragma once
#include "bsfwd.h"

#include "Solver/BSTargetInfo.h"
#include "Solver/SolverConfig.h"

namespace BSIPC
{
    namespace Energy
    {
        Float SeamValAt(const SeamPoint& seam, const SolverConfig& config, const BSSurface* surfI, const BSSurface* surfJ);

        Float SeamVal(const SolverConfig& config, const std::vector<BSTargetInfo>& infos);

        Mat<54, 1> SeamGradAt(const SeamPoint& seam, const SolverConfig& config,
            const BSSurface* surfI, const BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ);

        void SeamGrad(std::vector<BSTargetInfo> targets, const SolverConfig& config, BSIPC_OUT DMat& seamGrad);

        void SeamHessAt(const SeamPoint& seam, const SolverConfig& config,
            const BSSurface* surfI, const BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ,
            BSIPC_OUT DMat& localHess);

        void FillSeamGrad(const Mat<54, 1>& localGrad, BSIPC_OUT DMat& grad,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ, UInt offsetI, UInt offsetJ);

#if defined BSIPC_NUMERIC_TEST
        DMat NumericalSeamGradAt(const SeamPoint& seam, const SolverConfig& config,
            BSSurface* surfI, BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ);

        DMat NumericalSeamHessAt(const SeamPoint& seam, const SolverConfig& config,
            BSSurface* surfI, BSSurface* surfJ,
            std::vector<UInt> localIndicesI, std::vector<UInt> localIndicesJ);

        DMat NumericalSeamGrad(const SolverConfig& config, const std::vector<BSTargetInfo>& infos);
        DMat NumericalSeamHess(const SolverConfig& config, const std::vector<BSTargetInfo>& infos);
#endif
    }
}
