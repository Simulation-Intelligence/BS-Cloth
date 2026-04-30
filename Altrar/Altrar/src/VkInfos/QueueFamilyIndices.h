#pragma once
#include "atrfwd.h"

#include <array>

namespace ATR
{
    struct QueueFamilyIndices
    {
        static inline const size_t GRAPHICS = 0;        // graphic operations
        static inline const size_t PRESENT = 1;         // surface presentation support
        static inline const size_t TRANSFER = 2;        // transfer operations

        static inline const size_t COUNT = 3;

        std::array<std::optional<UInt>, QueueFamilyIndices::COUNT> indices;

        QueueFamilyIndices() 
        {
            for (auto& index : indices)
                index = std::nullopt;
        }

        Bool Complete()
        {
            for (auto& index : indices)
                if (index == std::nullopt)
                    return false;

            return true;
        }

        inline std::vector<UInt> UniqueIndices() 
        {
            std::vector<UInt> compactVec;
            for (size_t i = 0; i < QueueFamilyIndices::COUNT; i++)
                if (std::find(compactVec.begin(), compactVec.end(), (indices[i]).value()) == compactVec.end())
                    compactVec.push_back(indices[i].value());

            return compactVec;
        }

        inline Bool SeparateTransferQueue() { return indices[TRANSFER] != indices[GRAPHICS]; }
        inline Bool AllQueuesSame() { return indices[TRANSFER] == indices[GRAPHICS] && indices[GRAPHICS] == indices[PRESENT]; }

        friend std::ostream& operator<< (std::ostream& os, QueueFamilyIndices const& indices)
        {
            return os << 
                Format::item << "Graphics Queue Index: " << indices.indices[GRAPHICS].value() << "\n" <<
                Format::item << "Present Queue Index: " << indices.indices[PRESENT].value() << "\n" << 
                Format::item << "Transfer Queue Index: " << indices.indices[TRANSFER].value();
        }
    };
}