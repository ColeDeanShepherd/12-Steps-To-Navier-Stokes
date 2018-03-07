#include "Core.h"

int wrap(int value, int minValue, int maxValue)
{
    const auto numPossibleValues = (maxValue - minValue) + 1;

    while(value < minValue)
    {
        value += numPossibleValues;
    }

    while(value > maxValue)
    {
        value -= numPossibleValues;
    }

    return value;
}