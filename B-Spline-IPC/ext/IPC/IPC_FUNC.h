#pragma once
#ifndef _IPC_FUNC_H_
#define _IPC_FUNC_H_
//#include <vector>
//#include "collisionUtil.h"
//#include "Eigen/Eigen"
#include "fem3D.h"
#include "IPCdistanceFuncs.h"

namespace IPC
{
    //void Evaluate_GroundConstraintVals(const Ground& gd, const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
    //void Evaluate_SelfConstraintVals(const mesh3D& mesh, Eigen::VectorXd& vals, const int& offset);
    //double SelfConstraintVal(const mesh3D& mesh, const MMCVID& active);

    int solve_subIP(mesh3D& mesh, SpatialHash& sh, Ground& gd, double Kappa, float& time0, float& time1, float& time2, float& time3, float& time4, double& collisionNum);
    int IPC_Solver(int& stepId, mesh3D& mesh, SpatialHash& sh, Ground& gd);
    void ComputeXTilde(mesh3D& mesh);

    // Need to be called after properly filling in the mesh infos
    //  cf. FEMSimulator::buildModels()
    void BuildModel(mesh3D& mesh, SpatialHash& sh, Ground& gd);

    // Hardcode params here
    void InitSettings(mesh3D& mesh);

    void BuildFriction(mesh3D& mesh, const Ground& grd);

    void BuildCollisionSets(mesh3D& mesh, SpatialHash& sh, const Ground& gd, bool rehash);

    void SuggestKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass);
    void UpperBoundKappa(double& kappa, const double& Hhat, const double& bboxDiagSize2, const double& meanMass);
    void InitKappa(mesh3D& mesh, const Ground& grd, double& kappa);

    // WARNING the [moveDir] is the negative of the [searchDir] in BSIPC
    double EnvSelfCCDBound(mesh3D& mesh, SpatialHash& sh, Ground& gd, const std::vector<Eigen::Vector3d>& moveDir, double Kappa);

    // used in line search
    double PreLSIntersectionCCDBound(mesh3D& mesh, SpatialHash& sh, Ground& gd, double alpha,
        const std::vector<Eigen::Vector3d>& resultV0, const std::vector<Eigen::Vector3d>& searchDir, bool& intersected);

    // [preLSStepSize] is the step size output from [PreLSIntersectionCCDBound]
    double PostLSIntersectionCCDBound(double preLSStepSize, double stepSize, const std::vector<Eigen::Vector3d>& searchDir,
        mesh3D& mesh, SpatialHash& sh, Ground& gd, const std::vector<Eigen::Vector3d>& resultV0);

    void PostLineSearch(mesh3D& mesh, const Ground& grd, double alpha, double& kappa);

    // returns energy value
    double ComputeBarrierEnergyVal(const mesh3D& mesh, const Ground& gd, double Kappa);
    double ComputeFrictionEnergyVal(const mesh3D& mesh, const Ground& gd);

    // returns the number of collisions
    int ComputeBarrierFrictionGradHess(mesh3D& mesh, vector<Vector3d>& gradient, BHessian& BH, const Ground& grd);

    // Transform from DMat to vector<Vec3>
    Eigen::VectorXd TransformGrad(const std::vector<Eigen::Vector3d>& dir);
    std::vector<Eigen::Vector3d> TransformGrad(const Eigen::VectorXd& dir);

    void PreSolveAction(mesh3D& mesh);

    // If returns true, then break the loop, i.e. proceed to the next step
    bool PostSolveAction(Ground& gd, mesh3D& mesh);
}
#endif // !_IPC_FUNC_H_

