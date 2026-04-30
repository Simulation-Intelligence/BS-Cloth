#include "bspch.h"

#include "Inertia.h"

#include "Solver/Solver.h"

namespace BSIPC
{
    namespace Energy
    {
        Float InertiaVal(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos, const std::vector<Float>& massMatDiagEntries, Float thickness)
        {
            UInt size = static_cast<UInt>(estPos.size());

            Float inertiaVal = tbb::parallel_reduce(
                tbb::blocked_range<UInt>(0, size),
                0.0,
                [&](const tbb::blocked_range<UInt>& r, Float init) -> Float
                {
                    for (UInt i = r.begin(); i != r.end(); ++i)
                    {
                        Vec3 diff = estPos[i] - hypPos[i];
                        Float norm = diff.squaredNorm();
                        init += 0.5 * massMatDiagEntries[i] * norm;
                    }
                    return init;
                },
                std::plus<Float>()
            );

            return thickness * inertiaVal;
        }

        DMat InertiaGrad(const std::vector<Vec3>& hypPos, const std::vector<Vec3>& estPos, const std::vector<Float>& massMatDiagEntries, Float thickness)
        {
            DMat mat;
            UInt size = static_cast<UInt>(estPos.size());
            mat = DMat::Zero(3 * size, 1);

            tbb::parallel_for(
                0, static_cast<Int>(3 * size), 1, [&](Int idx)
                {
                    UInt i = idx / 3;
                    UInt j = idx % 3;

                    Vec3 diff = hypPos[i] - estPos[i];
                    Float curWeight = massMatDiagEntries[i];

                    mat(3 * i + j, 0) = curWeight * diff[j];
                }
            );

            return thickness * mat;
        }

        void FillInertiaHess(SpMatData& hessData, UInt offset, const std::vector<Float>& massMatDiagEntries, const SolverConfig& config, const Solver& solver)
        {
            UInt size = static_cast<UInt>(massMatDiagEntries.size());

            tbb::parallel_for(
                0, static_cast<Int>(3 * size), 1, [&](Int idx)
                {
                    UInt i = idx / 3;
                    UInt j = idx % 3;

                    UInt hessIndex = 3 * i + j;
                    if (solver.ShouldFixHessEntry(hessIndex, hessIndex))
                        hessData[offset + hessIndex] = SpMatEntry(hessIndex, hessIndex, 1);
                    else
                        hessData[offset + hessIndex] = SpMatEntry(hessIndex, hessIndex, config.thickness * massMatDiagEntries[i]);
                }
            );
        }
    }
}