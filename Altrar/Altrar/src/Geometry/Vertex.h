#pragma once

#include "atrfwd.h"

namespace ATR
{ 
    struct Vertex
    {
    public:
        Vertex(Vec3 pos, Vec3 normal) : pos(pos), normal(normal) {  }

        Vec3 pos;
        Vec3 normal;

        static VkVertexInputBindingDescription bindingDescription;
        static std::array<VkVertexInputAttributeDescription, 2> attributeDescriptions;

        inline String Info() const
        {
            String msg = "Pos: (" + std::to_string(pos.x) + ", " + std::to_string(pos.y) + ", " + std::to_string(pos.z) + 
                "), Normal: (" + std::to_string(normal.x) + ", " + std::to_string(normal.y) + ", " + std::to_string(normal.z) + ")";
            return msg;
        }
    };

    bool operator==(const Vertex& lhs, const Vertex& rhs);
}
