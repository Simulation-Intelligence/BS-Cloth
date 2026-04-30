#include "mIPC.h"

namespace IPC
{

    MMCVID::MMCVID(int a, int b, int c, int d)
    {
        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
    }
    MMCVID::MMCVID(int a)
    {
        data[0] = a;
        data[1] = data[2] = data[3] = -1;
    }

    void MMCVID::set(int a, int b, int c, int d)
    {
        data[0] = a;
        data[1] = b;
        data[2] = c;
        data[3] = d;
    }

}
