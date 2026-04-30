#pragma once

#include "bsfwd.h"

#include <chrono>
#include <unordered_map>

namespace BSIPC
{
    class Profiler
    {
    public:
        static void Start(const String& name);
        static void Check(const String& name);
        static void Update(const String& name);
        static void Stop(const String& name);
        static std::optional<UInt> GetTime(const String& name);

    private:
        inline static std::unordered_map<String, std::chrono::time_point<std::chrono::high_resolution_clock>> StartTimes;
        inline static std::unordered_map<String, UInt> Timings;           // in milliseconds
    };

#define BSIPC_PROFILE_SCOPE(name, scope) \
    Profiler::Start(name); \
    scope \
    Profiler::Check(name); \
    BSIPC_TRACE("{}: {}", name, Profiler::GetTime(name).value()); \
    Profiler::Stop(name);

#define BSIPC_PROFILE_START(name) \
    Profiler::Start(name);

#define BSIPC_PROFILE_END(name) \
    Profiler::Check(name); \
    BSIPC_TRACE("{}: {}", name, Profiler::GetTime(name).value()); \
    Profiler::Stop(name);
}
