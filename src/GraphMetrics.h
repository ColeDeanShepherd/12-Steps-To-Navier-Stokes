#pragma once

#include "Vector2d.h"

struct GraphMetrics
{
    double width, height;
    Vector2d pos;
    double minX, maxX;
    double minY, maxY;

    GraphMetrics() {}
};

Vector2d graphPointToPxPoint(const GraphMetrics& graphMetrics, const Vector2d& graphPoint);