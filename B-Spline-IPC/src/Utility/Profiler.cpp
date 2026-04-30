#include "bspch.h"
#include "Profiler.h"

namespace BSIPC
{
    void Profiler::Start(const String& name)
    {
        Profiler::StartTimes[name] = std::chrono::high_resolution_clock::now();
    }

    void Profiler::Check(const String& name)
    {
        if (Profiler::StartTimes.contains(name))
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto start = Profiler::StartTimes[name];  
            Profiler::Timings[name] = static_cast<UInt>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        }
    }

    void Profiler::Update(const String& name)
    {
        if (Profiler::StartTimes.contains(name))
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto start = Profiler::StartTimes[name];
            Profiler::Timings[name] = static_cast<UInt>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
            Profiler::StartTimes[name] = end;
        }
    }

    void Profiler::Stop(const String& name)
    {
        if (Profiler::StartTimes.contains(name))
        {
            auto end = std::chrono::high_resolution_clock::now();
            auto start = Profiler::StartTimes[name];
            Profiler::StartTimes.erase(name);
            Profiler::Timings[name] = static_cast<UInt>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
        }
    }

    std::optional<UInt> Profiler::GetTime(const String& name)
    {
        if (Profiler::Timings.contains(name))
            return Profiler::Timings[name];
        else
            return std::nullopt;
    }
}
