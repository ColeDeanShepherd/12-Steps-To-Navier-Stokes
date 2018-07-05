#pragma once

#include <string>
#include <vector>

template <typename T>
T clamp(const T value, const T min, const T max)
{
    if(value < min)
    {
        return min;
    }
    else if(value > max)
    {
        return max;
    }
    else
    {
        return value;
    }
}

int wrap(int value, const int minValue, const int maxValue);

std::vector<double> zeros(const size_t rowCount);
std::vector<std::vector<double>> zeros(const size_t rowCount, const size_t columnCount);
std::string toString(const std::vector<std::vector<double>>& matrix);