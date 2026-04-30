#include "atrpch.h"

#include "Loader/Config/Config.h"
#include "Renderer.h"

int main(int argc, char** argv)
{
    ATR::Config config;
    ATR::Renderer renderer(std::move(config));

    std::vector<ATR::Vertex> v = {
        {{-0.5f, -0.5f, 0.f}, {0.f, 0.f, 1.f}, {1.0f, 0.0f, 0.0f}},
        {{0.5f, -0.5f, 0.f}, {0.f, 0.f, 1.f}, {0.0f, 1.0f, 0.0f}},
        {{0.5f, 0.5f, 0.f}, {0.f, 0.f, 1.f}, {0.0f, 0.0f, 1.0f}},
        {{-0.5f, 0.5f, 0.f}, {0.f, 0.f, 1.f}, {1.0f, 1.0f, 1.0f}},

        {{-0.5f, -0.5f, -0.5f}, {0.f, 0.f, 1.f}, {1.0f, 0.0f, 0.0f}},
        {{0.5f, -0.5f, -0.5f}, {0.f, 0.f, 1.f}, {0.0f, 1.0f, 0.0f}},
        {{0.5f, 0.5f, -0.5f}, {0.f, 0.f, 1.f}, {0.0f, 0.0f, 1.0f}},
        {{-0.5f, 0.5f, -0.5f}, {0.f, 0.f, 1.f}, {1.0f, 1.0f, 1.0f}}
    };

    renderer.AddTriangle({v[0], v[1], v[2]});
    renderer.AddTriangle({v[2], v[3], v[0]});
    renderer.AddTriangle({v[4], v[5], v[6]});
    renderer.AddTriangle({v[6], v[7], v[4]});

    renderer.Run();

    return 0;
}
