#pragma once
#ifndef _IPC_TIME_STEP_FUNCS_
#define _IPC_TIME_STEP_FUNCS_
//#include <vector>
//#include "Eigen/Eigen"
#include "mesh.h"
#include "collisionUtil.h"
//void computeInjectiveStepSize_3d(const std::vector<std::vector<uint64_t>>& F,
//    const std::vector<Eigen::Vector3d>& x,
//    const std::vector<Eigen::Vector3d>& p,
//    double tol,
//    double slackness,
//    double* output);

namespace IPC
{

    void filterStepSize(const mesh3D& mesh, const std::vector<Eigen::Vector3d>& searchDir, double& stepSize);

    void Environment_largestFeasibleStepSize(const mesh3D& mesh, const Ground& gd,
        const std::vector<Eigen::Vector3d>& searchDir,
        double slackness,
        double& stepSize);

}

#endif // !_IPC_TIME_STEP_FUNCS_

