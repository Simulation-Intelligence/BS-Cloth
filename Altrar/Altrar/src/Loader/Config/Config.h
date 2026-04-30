#pragma once
#include "atrfwd.h"

namespace ATR
{
    struct Config
    {
        UInt width, height;
        Bool enableValidation;
        String location;
        std::vector<String> validationLayers;
        Bool saveScreenShots;
        UInt resumeAtFrame;

        Config();
        Config(const Config&) = default;
        Config(UInt, UInt, Bool, String, std::vector<String>, Bool, UInt);

        friend std::ostream& operator<< (std::ostream& os, Config const& config)
        {
            String layerStr = "";
            for (const auto& layer : config.validationLayers)
                layerStr += layer + "\n" + Format::subitem;
            layerStr = layerStr.substr(0, layerStr.size() - 4);

            return os <<
#if defined ATR_VERBOSE
                Format::item << "Verbose Mode on.\n" <<
#else
                Format::item << "Verbose Mode off.\n" <<
#endif
                Format::item << "Width: " << config.width << ", " << "Height: " << config.height << "\n" <<
                Format::item << "Enable Validation: " << config.enableValidation << "\n" <<
                Format::item << "Validation Layers: \n" << Format::subitem << layerStr;
        }
        
        static Config defaultConfig;
    };

}