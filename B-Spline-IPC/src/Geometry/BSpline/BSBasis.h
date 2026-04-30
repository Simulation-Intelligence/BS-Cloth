#pragma once

#pragma warning(disable: 4702)

#include "bsfwd.h"

#include "Utility/Math/MathHelpers.h"

namespace BSIPC
{
    using BSBasisFn = std::function<Float(Float)>;

    using BSEvalFn1 = std::function<Float(Float)>;
    using BSEvalFn2 = std::function<Float(Float, Float)>;

    namespace BSBasis
    {
        inline Float LeftBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 1, "Left Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return (x - 1.) * (x - 1.);

            BSIPC_WARN("B-Spline Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float LeftInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 2, "Left Boundary B-Spline derivative undefined for x > 2: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return -3. / 2. * x * x + 2. * x;
            else if (InClosedInterval(x, 1., 2.))
                return 0.5f * x * x - 2. * x + 2.;

            BSIPC_WARN("B-Spline Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float Interior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 3, "Interior B-Spline derivative undefined for x > 3: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return 0.5f * x * x;
            else if (InClosedInterval(x, 1., 2.))
                return -x * x + 3. * x - 3. / 2.;
            else if (InClosedInterval(x, 2., 3.))
                return 0.5f * (x - 3.) * (x - 3.);

            BSIPC_WARN("B-Spline Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 2, "Right Interior B-Spline derivative undefined for x > 2: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return 0.5f * x * x;
            else if (InClosedInterval(x, 1., 2.))
                return -3. / 2. * x * x + 4. * x - 2.;

            BSIPC_WARN("B-Spline Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 1, "Right Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return x * x;

            BSIPC_WARN("B-Spline Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        const std::array<BSBasisFn, 5> BasisFns = { LeftBoundary, LeftInterior, Interior, RightInterior, RightBoundary };
    }

    namespace BSBasisDerivative
    {
        inline Float LeftBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 1, "Left Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return 2. * x - 2.;
            else
                return 0.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float LeftInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 2, "Left Boundary B-Spline derivative undefined for x > 2: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return -3. * x + 2.;
            else if (InClosedInterval(x, 1., 2.))
                return x - 2.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float Interior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 3, "Interior B-Spline derivative undefined for x > 3: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return x;
            else if (InClosedInterval(x, 1., 2.))
                return -2. * x + 3.;
            else if (InClosedInterval(x, 2., 3.))
                return x - 3.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 2, "Right Interior B-Spline derivative undefined for x > 2: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return x;
            else if (InClosedInterval(x, 1., 2.))
                return -3. * x + 4.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 1, "Right Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return 2. * x;
            
            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        const std::array<BSBasisFn, 5> BasisFns = { LeftBoundary, LeftInterior, Interior, RightInterior, RightBoundary };
    }

    namespace BSBasisSecondDerivative
    {
        inline Float LeftBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
            BSIPC_ASSERT(x <= 1, "Left Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

            if (InClosedInterval(x, 0., 1.))
                return 2.;
            else
                return 0.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float LeftInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Left Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
                BSIPC_ASSERT(x <= 2, "Left Boundary B-Spline derivative undefined for x > 2: " + std::to_string(x))

                if (InClosedInterval(x, 0., 1.))
                    return -3.;
                else if (InClosedInterval(x, 1., 2.))
                    return 1.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float Interior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
                BSIPC_ASSERT(x <= 3, "Interior B-Spline derivative undefined for x > 3: " + std::to_string(x))

                if (InClosedInterval(x, 0., 1.))
                    return 1.;
                else if (InClosedInterval(x, 1., 2.))
                    return -2.;
                else if (InClosedInterval(x, 2., 3.))
                    return 1.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightInterior(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Interior B-Spline derivative undefined for negative x: " + std::to_string(x))
                BSIPC_ASSERT(x <= 2, "Right Interior B-Spline derivative undefined for x > 2: " + std::to_string(x))

                if (InClosedInterval(x, 0., 1.))
                    return 1.;
                else if (InClosedInterval(x, 1., 2.))
                    return -3.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        inline Float RightBoundary(Float x)
        {
            BSIPC_ASSERT(x >= 0, "Right Boundary B-Spline derivative undefined for negative x: " + std::to_string(x))
                BSIPC_ASSERT(x <= 1, "Right Boundary B-Spline derivative undefined for x > 1: " + std::to_string(x))

                if (InClosedInterval(x, 0., 1.))
                    return 2.;

            BSIPC_WARN("B-Spline Derivative Evaluation Falling Through; X: {}", x);
            return std::numeric_limits<Float>::max();
        }

        const std::array<BSBasisFn, 5> BasisFns = { LeftBoundary, LeftInterior, Interior, RightInterior, RightBoundary };
    }
}

#pragma warning(default: 4702)
