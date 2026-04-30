#pragma once

#include "bsfwd.h"

namespace BSIPC
{
    enum TransformationType
    {
        Translation,
        Rotation,
        Scaling
    };

    /** This is to emulate a union (additive-enum) type.
     *  
     *  type = Translation:     content is the translational vector
     *  type = Rotation:        content is the scaling in degree (applied in the sequence of x/y/z)
     *  type = Scaling:         content is the scaling coefficient along each axis
     * 
     *  For rotation, currently supporting only rotation around x/y/z axis.
     */
    struct LinearTransformation
    {
        TransformationType type;
        Vec3 content;

        inline void ApplyTo(BSIPC_OUT Vec3& vec) const;
    };

    inline void LinearTransformation::ApplyTo(BSIPC_OUT Vec3& vec) const
    {
        Vec3 content = this->content;

        if (this->type == TransformationType::Translation)
            vec += content;
        else if (this->type == TransformationType::Rotation)
        {
            Vec3 contentInRad = content * Const::PI / 180;
            Mat<3, 3> rotateX, rotateY, rotateZ;
            rotateX << 1, 0, 0,
                0, cos(contentInRad[0]), -sin(contentInRad[0]),
                0, sin(contentInRad[0]), cos(contentInRad[0]);
            rotateY << cos(contentInRad[1]), 0, sin(contentInRad[1]),
                0, 1, 0,
                -sin(contentInRad[1]), 0, cos(contentInRad[1]);
            rotateZ << cos(contentInRad[2]), -sin(contentInRad[2]), 0,
                sin(contentInRad[2]), cos(contentInRad[2]), 0,
                0, 0, 1;

            Mat<3, 3> rotate = rotateZ * rotateY * rotateX;
            vec = rotate * vec;
        }
        else if (this->type == TransformationType::Scaling)
        {
            vec = Vec3(
                vec[0] * content[0],
                vec[1] * content[1],
                vec[2] * content[2]
            );
        }
        else
        {
            BSIPC_WARN("Transformation Type Fallthrough.");
        }
    }
}