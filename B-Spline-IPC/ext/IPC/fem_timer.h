#pragma once
#ifndef FEM_TIMER_H
#define FEM_TIMER_H
// 
#include <ctime>

namespace IPC
{

    class HighResolutionTimer
    {
    public:
        virtual void set_start() = 0;
        virtual void set_end() = 0;
        virtual float get_millisecond() = 0;
    };

#ifdef WIN32
#include <windows.h>
    class HighResolutionTimerForWin : public HighResolutionTimer
    {
    public:

        HighResolutionTimerForWin() {
            QueryPerformanceFrequency(&freq_);
            start_.QuadPart = 0;
            end_.QuadPart = 0;
        }

        void set_start() {
            QueryPerformanceCounter(&start_);
        }

        void set_end() {
            QueryPerformanceCounter(&end_);
        }

        float get_millisecond() {
            return static_cast<float>((end_.QuadPart - start_.QuadPart) * 1000 / (float)freq_.QuadPart);
        }

    private:
        LARGE_INTEGER freq_;
        LARGE_INTEGER start_, end_;
    };
#else
    class HighResolutionTimerForWin : public HighResolutionTimer
    {
    public:

        HighResolutionTimerForWin() {
            start_ = 0;
            end_ = 0;
        }

        void set_start() {
            struct timespec t;
            std::timespec_get(&t, TIME_UTC);
            start_ = t.tv_sec * 1e3 + t.tv_nsec * 1e-6;
        }

        void set_end() {
            struct timespec t;
            std::timespec_get(&t, TIME_UTC);
            end_ = t.tv_sec * 1e3 + t.tv_nsec * 1e-6;
        }

        float get_millisecond() {
            // return static_cast<float>((end_.QuadPart - start_.QuadPart) * 1000 / (float)freq_.QuadPart);
            return static_cast<float>(end_ - start_);   // ms
        }

    private:
        double start_, end_;
    };
#endif

}
#endif
