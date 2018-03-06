#pragma once

namespace std
{
    template <typename T>
    T clamp(T value, T min, T max)
    {
        if(value < min)
        {
            return min;
        } else if(value > max)
        {
            return max;
        } else
        {
            return value;
        }
    }
}

const int WINDOW_WIDTH = 640;
const int WINDOW_HEIGHT = 480;