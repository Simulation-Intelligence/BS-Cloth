#include "bspch.h"
#include "BSSurfaceMgr.h"

#include "Solver.h"

#include "Utility/EigenHelper.h"

namespace BSIPC
{
    void BSSurfaceMgr::ParseFromConfig(const SolverConfig& config)
    {
        UInt index = 0;
        for (const BSSurfaceInfo& info : config.surfacesInfo)
        {
            BSSurface surf(config, index, info);
            this->surfaces.push_back(surf);
            if (info.initVelocities.empty())
                this->initVels.push_back(info.initVel);
            else
                this->initVels.push_back(info.initVelocities);
            ++index;
        }
    }

    void BSSurfaceMgr::LinkSolver(Solver& solver)
    {
        for (UInt i = 0; i != this->surfaces.size(); ++i)
        {
            if (std::holds_alternative<Vec3>(this->initVels[i]))
                solver.AddTarget(&this->surfaces[i], std::get<Vec3>(this->initVels[i]));
            else if (std::holds_alternative<std::vector<Vec3>>(this->initVels[i]))
            {
                const std::vector<Vec3>& initVelsVec = std::get<std::vector<Vec3>>(this->initVels[i]);
                solver.AddTarget(&this->surfaces[i], initVelsVec);
            }
            else
                BSIPC_CRITICAL("BSSurfaceMgr::LinkSolver variant enumeration fallthrough!");
            this->surfaces[i].BindSolver(&solver);
        }
    }
}
