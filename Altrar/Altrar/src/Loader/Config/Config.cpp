#include "atrpch.h"

#include "Config.h"

namespace ATR
{
    //#define LOAD_DATA_FROM_YAML(var, node, tag, TYPE) \
    //    if (node.contains(#tag)) \
    //        var = node[#tag].get_value<TYPE>(); \
    //    else \
    //        throw fkyaml::exception(("missing tag " + std::string(#tag)).c_str());
    //
    //#define LOAD_DATA_FROM_YAML_NOERROR(var, node, tag, TYPE) \
    //    if (node.contains(#tag)) \
    //        var = root[#tag].get_value<TYPE>(); \
    //    else \
    //        ATR_PRINT("[INFO] tag \"" + String(#tag) + "\" not specified.");
    //
    //#define LOAD_DEF_DATA_FROM_YAML(var, node, tag, TYPE) \
    //    TYPE var; \
    //    LOAD_DATA_FROM_YAML(var, node, tag, TYPE)
    //
    //#define LOAD_DEF_DATA_FROM_YAML_NOERROR(var, node, tag, TYPE) \
    //    TYPE var; \
    //    LOAD_DATA_FROM_YAML_NOERROR(var, node, tag, TYPE)
    //
    //#define LOAD_NODE_FROM_YAML(node, root, tag) \
    //    auto node = root; \
    //    if (root.contains(#tag)) \
    //        node = root[#tag]; \
    //    else \
    //        throw fkyaml::exception(("missing tag " + std::string(#tag)).c_str());
    //
    //#define LOAD_NODE_FROM_YAML_NOERROR(node, root, tag) \
    //    auto node = root; \
    //    if (root.contains(#tag)) \
    //        node = root[#tag]; \
    //    else \
    //        ATR_PRINT("[INFO] node \"" + String(#tag) + "\" not specified.");

    Config::Config() :
        width(800), height(600),
        enableValidation(true),
        location(""),
        validationLayers({ "VK_LAYER_KHRONOS_validation" }),
        saveScreenShots(false),
        resumeAtFrame(0)
    {   }

    Config::Config
    (
        UInt width, UInt Height,
        Bool validation,
        String path,
        std::vector<String> validationLayers,
        Bool saveScreenShots,
        UInt resumeAtFrame = 0
    ) :
        width(width), height(Height),
        enableValidation(validation),
        location(path),
        validationLayers(validationLayers),
        saveScreenShots(saveScreenShots),
        resumeAtFrame(resumeAtFrame)
    {  }

    Config Config::defaultConfig
    (
        800, 800,
        true,
        "..\\Altrar\\Altrar",
        { "VK_LAYER_KHRONOS_validation" },
        false,
        0
    );
}