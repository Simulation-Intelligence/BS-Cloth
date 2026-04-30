#pragma once
#include "Eigen/Eigen"

namespace IPC
{
    using namespace Eigen;
    double point_triangle_ccd(
        const Vector3d& _p,
        const Vector3d& _t0,
        const Vector3d& _t1,
        const Vector3d& _t2,
        const Vector3d& _dp,
        const Vector3d& _dt0,
        const Vector3d& _dt1,
        const Vector3d& _dt2,
        const unsigned int _pidx, const unsigned int _t0idx, const unsigned int _t1idx, const unsigned int _t2idx,
        double eta, double thickness);

    double edge_edge_ccd(
        const Vector3d& _ea0,
        const Vector3d& _ea1,
        const Vector3d& _eb0,
        const Vector3d& _eb1,
        const Vector3d& _dea0,
        const Vector3d& _dea1,
        const Vector3d& _deb0,
        const Vector3d& _deb1,
        double eta, double thickness);

    bool point_triangle_ccd_broadphase(
        const Vector3d& p,
        const Vector3d& t0,
        const Vector3d& t1,
        const Vector3d& t2,
        const Vector3d& dp,
        const Vector3d& dt0,
        const Vector3d& dt1,
        const Vector3d& dt2,
        double dist);

    bool edge_edge_ccd_broadphase(
        const Vector3d& ea0,
        const Vector3d& ea1,
        const Vector3d& eb0,
        const Vector3d& eb1,
        const Vector3d& dea0,
        const Vector3d& dea1,
        const Vector3d& deb0,
        const Vector3d& deb1,
        double dist);
}
