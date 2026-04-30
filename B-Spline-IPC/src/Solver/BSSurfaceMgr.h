#pragma once

#include "bsfwd.h"
#include <variant>

#include "Geometry/BSpline/BSSurface.h"

#include "Solver/SolverConfig.h"

namespace BSIPC
{
    class BSSurfaceMgr
    {
    public:
        void ParseFromConfig(const SolverConfig& config);
        void LinkSolver(BSIPC_OUT Solver& solver);

    private:
        std::vector<BSSurface> surfaces;
        std::vector<std::variant<Vec3, std::vector<Vec3>>> initVels;
    };
}
