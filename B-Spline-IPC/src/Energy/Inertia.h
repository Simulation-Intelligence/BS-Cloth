#pragma once
#include "bsfwd.h"

#include "Solver/SolverConfig.h"
#include "Solver/BSTargetInfo.h"

namespace BSIPC
{
    class Solver;

    namespace Energy
    {
        /* ***************************** Inertia ***************************** */

        Float InertiaVal(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos, const std::vector<Float>& massMatDiagEntries, Float thickness);

        // Should return a 3n-by-1 column vector
        DMat InertiaGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos, const std::vector<Float>& massMatDiagEntries, Float thickness);

        // Addition on SpMat is expensive. Writes to the global [hessData] directly, starting from index [offset]
        void FillInertiaHess(SpMatData& hessData, UInt offset, const std::vector<Float>& massMatDiagEntries, const SolverConfig& config, const Solver& solver);
    }
}