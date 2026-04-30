#pragma once
#include "bsfwd.h"

#include <fstream>

#include "Core/Log.h"
#include "fkyaml/node.hpp"

namespace BSIPC
{
    template <typename T>
    inline void LoadDataFromYaml(T& var, const fkyaml::node& node, const char* tag)
    {
        if (node.contains(tag))
            var = node[tag].get_value<T>();
        else
            throw fkyaml::exception(("missing tag " + std::string(tag)).c_str());
    }

    template <typename T>
    inline bool LoadDataFromYamlOptional(T& var, const fkyaml::node& node, const char* tag)
    {
        if (node.contains(tag))
        {
            var = node[tag].get_value<T>();
            return true;
        }
        else
        {
            BSIPC_INFO("tag \"" + String(tag) + "\" not specified.");
            return false;
        }
    }
   
    inline void LoadNodeFromYaml(fkyaml::node& node, const fkyaml::node& parent, const char* tag)
    {
        if (parent.contains(tag))
            node = parent[tag];
        else
            throw fkyaml::exception(("missing tag " + std::string(tag)).c_str());
    }

    inline bool LoadNodeFromYamlOptional(fkyaml::node& node, const fkyaml::node& parent, const char* tag)
    {
        if (parent.contains(tag))
        {
            node = parent[tag];
            return true;
        }
        else
        {
            BSIPC_INFO("node \"" + String(tag) + "\" not specified.");
            return false;
        }
    }

    inline String ToStr(const fkyaml::node& node)
    {
        std::stringstream ss;
        ss << node;
        return ss.str();
    }
}
