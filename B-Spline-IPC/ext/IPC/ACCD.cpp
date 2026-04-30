#include "ACCD.h"

namespace IPC 
{
    using namespace Eigen;

    int point_triangle_distance_type(
        const Vector3d& p,
        const Vector3d& t0,
        const Vector3d& t1,
        const Vector3d& t2)
    {
        Eigen::Matrix<double, 2, 3> basis;
        basis.row(0) = (t1 - t0).transpose();
        basis.row(1) = (t2 - t0).transpose();

        const Eigen::Matrix<double, 3, 1> nVec = basis.row(0).cross(basis.row(1));

        Eigen::Matrix<double, 2, 3> param;

        basis.row(1) = basis.row(0).cross(nVec);
        param.col(0) = (basis * basis.transpose()).ldlt().solve(basis * (p - t0));
        if (param(0, 0) > 0.0 && param(0, 0) < 1.0 && param(1, 0) >= 0.0) {
            return 3; // PE t0t1
        }
        else {
            basis.row(0) = (t2 - t1).transpose();

            basis.row(1) = basis.row(0).cross(nVec);
            param.col(1) = (basis * basis.transpose()).ldlt().solve(basis * (p - t1));
            if (param(0, 1) > 0.0 && param(0, 1) < 1.0 && param(1, 1) >= 0.0) {
                return 4; // PE t1t2
            }
            else {
                basis.row(0) = (t0 - t2).transpose();

                basis.row(1) = basis.row(0).cross(nVec);
                param.col(2) = (basis * basis.transpose()).ldlt().solve(basis * (p - t2));
                if (param(0, 2) > 0.0 && param(0, 2) < 1.0 && param(1, 2) >= 0.0) {
                    return 5; // PE t2t0
                }
                else {
                    if (param(0, 0) <= 0.0 && param(0, 2) >= 1.0) {
                        return 0; // PP t0
                    }
                    else if (param(0, 1) <= 0.0 && param(0, 0) >= 1.0) {
                        return 1; // PP t1
                    }
                    else if (param(0, 2) <= 0.0 && param(0, 1) >= 1.0) {
                        return 2; // PP t2
                    }
                    else {
                        return 6; // PT
                    }
                }
            }
        }
    }

    int edge_edge_distance_type(
        const Vector3d& ea0,
        const Vector3d& ea1,
        const Vector3d& eb0,
        const Vector3d& eb1)
    {
        auto u = ea1 - ea0;
        auto v = eb1 - eb0;
        auto w = ea0 - eb0;
        double a = u.squaredNorm(); // always >= 0
        double b = u.dot(v);
        double c = v.squaredNorm(); // always >= 0
        double d = u.dot(w);
        double e = v.dot(w);
        double D = a * c - b * b; // always >= 0
        double tD = D; // tc = tN / tD, default tD = D >= 0
        double sN, tN;

        int defaultCase = 8;

        // compute the line parameters of the two closest points
        sN = (b * e - c * d);
        if (sN <= 0.0) { // sc < 0 => the s=0 edge is visible
            tN = e;
            tD = c;
            defaultCase = 2;
        }
        else if (sN >= D) { // sc > 1  => the s=1 edge is visible
            tN = e + b;
            tD = c;
            defaultCase = 5;
        }
        else {
            tN = (a * e - b * d);
            if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() < 1.0e-20 * a * c)) {
                // if (tN > 0.0 && tN < tD && (u.cross(v).dot(w) == 0.0 || u.cross(v).squaredNorm() == 0.0)) {
                // std::cout << u.cross(v).squaredNorm() / (a * c) << ": " << sN << " " << D << ", " << tN << " " << tD << std::endl;
                // avoid coplanar or nearly parallel EE
                if (sN < D / 2) {
                    tN = e;
                    tD = c;
                    defaultCase = 2;
                }
                else {
                    tN = e + b;
                    tD = c;
                    defaultCase = 5;
                }
            }
            // else defaultCase stays as 8
        }

        if (tN <= 0.0) { // tc < 0 => the t=0 edge is visible
            // recompute sc for this edge
            if (-d <= 0.0) {
                return 0;
            }
            else if (-d >= a) {
                return 3;
            }
            else {
                return 6;
            }
        }
        else if (tN >= tD) { // tc > 1  => the t=1 edge is visible
            // recompute sc for this edge
            if ((-d + b) <= 0.0) {
                return 1;
            }
            else if ((-d + b) >= a) {
                return 4;
            }
            else {
                return 7;
            }
        }

        return defaultCase;
    }

    double point_point_distance(const Vector3d& a, const Vector3d& b)
    {
        return (a - b).squaredNorm();
    }

    double point_edge_distance(const Vector3d& p, const Vector3d& e0, const Vector3d& e1)
    {
        return (e0 - p).cross(e1 - p).squaredNorm() / (e1 - e0).squaredNorm();
    }

    double point_triangle_distance(const Vector3d& p, const Vector3d& t0, const Vector3d& t1, const Vector3d& t2)
    {
        const Eigen::Matrix<double, 3, 1> b = (t1 - t0).cross(t2 - t0);
        double aTb = (p - t0).dot(b);
        return aTb * aTb / b.squaredNorm();
    }

    double edge_edge_distance(const Vector3d& ea0, const Vector3d& ea1, const Vector3d& eb0, const Vector3d& eb1)
    {
        auto b = (ea1 - ea0).cross(eb1 - eb0);
        double aTb = (eb0 - ea0).dot(b);
        return aTb * aTb / b.squaredNorm();
    }

    bool edge_edge_ccd_broadphase(
        const Vector3d& ea0,
        const Vector3d& ea1,
        const Vector3d& eb0,
        const Vector3d& eb1,
        const Vector3d& dea0,
        const Vector3d& dea1,
        const Vector3d& deb0,
        const Vector3d& deb1,
        double dist)
    {
        auto max_a = ea0.array().max(ea1.array()).max((dea0).array()).max((dea1).array());
        auto min_a = ea0.array().min(ea1.array()).min((dea0).array()).min((dea1).array());
        auto max_b = eb0.array().max(eb1.array()).max((deb0).array()).max((deb1).array());
        auto min_b = eb0.array().min(eb1.array()).min((deb0).array()).min((deb1).array());
        if ((min_a - max_b > dist).any() || (min_b - max_a > dist).any())
            return false;
        else
            return true;
    }

    bool point_triangle_ccd_broadphase(
        const Vector3d& p,
        const Vector3d& t0,
        const Vector3d& t1,
        const Vector3d& t2,
        const Vector3d& dp,
        const Vector3d& dt0,
        const Vector3d& dt1,
        const Vector3d& dt2,
        double dist)
    {
        auto max_p = p.array().max((dp).array());
        auto min_p = p.array().min((dp).array());
        auto max_tri = t0.array().max(t1.array()).max(t2.array()).max((dt0).array()).max((dt1).array()).max((dt2).array());
        auto min_tri = t0.array().min(t1.array()).min(t2.array()).min((dt0).array()).min((dt1).array()).min((dt2).array());
        if ((min_p - max_tri > dist).any() || (min_tri - max_p > dist).any())
            return false;
        else
            return true;
    }

    double edge_edge_distance_unclassified(
        const Vector3d& ea0,
        const Vector3d& ea1,
        const Vector3d& eb0,
        const Vector3d& eb1)
    {
        switch (edge_edge_distance_type(ea0, ea1, eb0, eb1)) {
        case 0:
            return point_point_distance(ea0, eb0);
        case 1:
            return point_point_distance(ea0, eb1);
        case 2:
            return point_edge_distance(ea0, eb0, eb1);
        case 3:
            return point_point_distance(ea1, eb0);
        case 4:
            return point_point_distance(ea1, eb1);
        case 5:
            return point_edge_distance(ea1, eb0, eb1);
        case 6:
            return point_edge_distance(eb0, ea0, ea1);
        case 7:
            return point_edge_distance(eb1, ea0, ea1);
        case 8:
            return edge_edge_distance(ea0, ea1, eb0, eb1);
        default:
            return std::numeric_limits<double>::max();
        }
    }

    double point_triangle_distance_unclassified(
        const Vector3d& p,
        const Vector3d& t0,
        const Vector3d& t1,
        const Vector3d& t2)
    {
        switch (point_triangle_distance_type(p, t0, t1, t2)) {
        case 0:
            return point_point_distance(p, t0);
        case 1:
            return point_point_distance(p, t1);
        case 2:
            return point_point_distance(p, t2);
        case 3:
            return point_edge_distance(p, t0, t1);
        case 4:
            return point_edge_distance(p, t1, t2);
        case 5:
            return point_edge_distance(p, t2, t0);
        case 6:
            return point_triangle_distance(p, t0, t1, t2);
        default:
            return std::numeric_limits<double>::max();
        }
    }

    double edge_edge_ccd(
        const Vector3d& _ea0,
        const Vector3d& _ea1,
        const Vector3d& _eb0,
        const Vector3d& _eb1,
        const Vector3d& _dea0,
        const Vector3d& _dea1,
        const Vector3d& _deb0,
        const Vector3d& _deb1,
        double eta, double thickness)
    {
        Vector3d ea0 = _ea0, ea1 = _ea1, eb0 = _eb0, eb1 = _eb1, dea0 = _dea0, dea1 = _dea1, deb0 = _deb0, deb1 = _deb1;
        Vector3d mov = (dea0 + dea1 + deb0 + deb1) / 4;
        dea0 -= mov;
        dea1 -= mov;
        deb0 -= mov;
        deb1 -= mov;
        double max_disp_mag = std::sqrt(std::max(dea0.squaredNorm(), dea1.squaredNorm())) + std::sqrt(std::max(deb0.squaredNorm(), deb1.squaredNorm()));
        if (max_disp_mag == 0)
            return 1.0;

        double dist2_cur = edge_edge_distance_unclassified(ea0, ea1, eb0, eb1);
        //printf("dist2_cur = %f\n", dist2_cur);
        double dFunc = dist2_cur - thickness * thickness;
        if (dFunc <= 0) {
            std::vector<double> dists{ (ea0 - eb0).squaredNorm(), (ea0 - eb1).squaredNorm(), (ea1 - eb0).squaredNorm(), (ea1 - eb1).squaredNorm() };
            dist2_cur = *std::min_element(dists.begin(), dists.end());
            dFunc = dist2_cur - thickness * thickness;
        }
        double dist_cur = std::sqrt(dist2_cur);
        double gap = eta * dFunc / (dist_cur + thickness);
        double toc = 0.0;
        int iter_cnt = 0;

        // Guards against infinite loops
        const int MAX_ITERATIONS = 10000;
        const double MIN_PROGRESS = 1e-10;
        const double TOC_EPSILON = 1e-8;

        while (true) {
            // Guard 1: Maximum iteration limit
            if (iter_cnt >= MAX_ITERATIONS) {
                printf("Warning: edge-edge CCD max iterations reached\n");
                return 1.0;
            }

            double toc_lower_bound = (1 - eta) * dFunc / ((dist_cur + thickness) * max_disp_mag);

            // Guard 2: Check for numerical issues or lack of progress
            if (toc_lower_bound < MIN_PROGRESS || !std::isfinite(toc_lower_bound)) {
                // fprintf(stderr, "Warning: edge-edge CCD numerical issue detected\n");
                return 1.0; // Return current time or 1.0 depending on desired behavior
            }

            ea0 += toc_lower_bound * dea0;
            ea1 += toc_lower_bound * dea1;
            eb0 += toc_lower_bound * deb0;
            eb1 += toc_lower_bound * deb1;
            dist2_cur = edge_edge_distance_unclassified(ea0, ea1, eb0, eb1);
            dFunc = dist2_cur - thickness * thickness;
            if (dFunc <= 0) {
                std::vector<double> dists{ (ea0 - eb0).squaredNorm(), (ea0 - eb1).squaredNorm(), (ea1 - eb0).squaredNorm(), (ea1 - eb1).squaredNorm() };
                dist2_cur = *std::min_element(dists.begin(), dists.end());
                dFunc = dist2_cur - thickness * thickness;
            }
            dist_cur = std::sqrt(dist2_cur);

            if (toc > TOC_EPSILON && (dFunc / (dist_cur + thickness) < gap)) {
                break;
            }

            toc += toc_lower_bound;
            if (toc >= 1.0 - TOC_EPSILON)
                return 1.0;

            ++iter_cnt;
        }
        return toc;
    }

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
        double eta, double thickness)
    {
        Vector3d p = _p, t0 = _t0, t1 = _t1, t2 = _t2, dp = _dp, dt0 = _dt0, dt1 = _dt1, dt2 = _dt2;
        Vector3d mov = (dt0 + dt1 + dt2 + dp) / 4;
        dt0 -= mov;
        dt1 -= mov;
        dt2 -= mov;
        dp -= mov;
        std::vector<double> disp_mag2_vec{ dt0.squaredNorm(), dt1.squaredNorm(), dt2.squaredNorm() };
        double max_disp_mag = dp.norm() + std::sqrt(*std::max_element(disp_mag2_vec.begin(), disp_mag2_vec.end()));
        if (max_disp_mag == 0)
            return 1.0;

        double dist2_cur = point_triangle_distance_unclassified(p, t0, t1, t2);
        double dist_cur = std::sqrt(dist2_cur);
        double gap = eta * (dist2_cur - thickness * thickness) / (dist_cur + thickness);
        double toc = 0.0;
        int iter_cnt = 0;

        const int MAX_ITERATIONS = 10000;
        const double MIN_PROGRESS = 1e-10;
        const double TOC_EPSILON = 1e-8;

        while (true) {
            if (iter_cnt >= MAX_ITERATIONS) {
                printf("Warning: point-triangle CCD max iterations reached || %u, (%u, %u, %u)\n", _pidx, _t0idx, _t1idx, _t2idx);
                return 1.0;
            }

            double toc_lower_bound = (1 - eta) * (dist2_cur - thickness * thickness)
                / ((dist_cur + thickness) * max_disp_mag);

            if (toc_lower_bound < MIN_PROGRESS || !std::isfinite(toc_lower_bound)) {
                return 1.0; // Numerical issue - return current time
            }

            p += toc_lower_bound * dp;
            t0 += toc_lower_bound * dt0;
            t1 += toc_lower_bound * dt1;
            t2 += toc_lower_bound * dt2;

            dist2_cur = point_triangle_distance_unclassified(p, t0, t1, t2);
            dist_cur = std::sqrt(dist2_cur);

            if (toc > TOC_EPSILON && ((dist2_cur - thickness * thickness) / (dist_cur + thickness) < gap)) {
                break;
            }

            toc += toc_lower_bound;
            if (toc >= 1.0 - TOC_EPSILON) {
                return 1.0;
            }
            ++iter_cnt;
        }
        return toc;
    }

}