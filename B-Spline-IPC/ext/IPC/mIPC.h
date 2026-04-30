#pragma once
#ifndef _M_IPC_H_
#define _M_IPC_H_
#include <array>

namespace IPC
{

class MMCVID {
public:
    std::array<int, 4> data;

public:
    MMCVID(int a, int b, int c, int d);
    MMCVID(int a = -1);

    void set(int a, int b, int c, int d);

    int& operator[](int i)
    {
        return data[i];
    }

    const int& operator[](int i) const
    {
        return data[i];
    }

    friend bool operator==(const MMCVID& lhs, const MMCVID& rhs)
    {
        if ((lhs[0] == rhs[0]) && (lhs[1] == rhs[1]) && (lhs[2] == rhs[2]) && (lhs[3] == rhs[3])) {
            return true;
        }
        else {
            return false;
        }
    }

    friend bool operator<(const MMCVID& lhs, const MMCVID& rhs)
    {
        if (lhs[0] < rhs[0]) {
            return true;
        }
        else if (lhs[0] == rhs[0]) {
            if (lhs[1] < rhs[1]) {
                return true;
            }
            else if (lhs[1] == rhs[1]) {
                if (lhs[2] < rhs[2]) {
                    return true;
                }
                else if (lhs[2] == rhs[2]) {
                    if (lhs[3] < rhs[3]) {
                        return true;
                    }
                    else {
                        return false;
                    }
                }
                else {
                    return false;
                }
            }
            else {
                return false;
            }
        }
        else {
            return false;
        }
    }
};

}

#endif // !_M_IPC_H_

