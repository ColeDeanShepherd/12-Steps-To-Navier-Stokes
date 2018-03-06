#pragma once

#include <vector>
#include "Vector2d.h"

double gradient1stOrderBackwardDiff(const std::vector<double>& f, const size_t i, const double dx)
{
    return (f[i] - f[i - 1]) / dx;
}
double gradient2ndOrderCentralDiff(const std::vector<double>& f, const size_t i, const double dx)
{
    return (f[i + 1] - (2 * f[i]) + f[i - 1]) / (dx * dx);
}
Vector2d gradient1stOrderCentralDiff(
    const std::vector<std::vector<double>>& f,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto dfdx = (f[xIndex + 1][yIndex] - f[xIndex - 1][yIndex]) / (2 * dx);
    const auto dfdy = (f[xIndex][yIndex + 1] - f[xIndex][yIndex - 1]) / (2 * dy);

    return Vector2d(dfdx, dfdy);
}
double divergence1stOrderBackwardDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto dvxdx = (v[xIndex][yIndex].x - v[xIndex - 1][yIndex].x) / dx;
    const auto dvydy = (v[xIndex][yIndex].y - v[xIndex][yIndex - 1].y) / dy;

    return dvxdx + dvydy;
}
Vector2d laplacian2ndOrderCentralDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto d2vxdx = (v[xIndex + 1][yIndex].x - (2 * v[xIndex][yIndex].x) + v[xIndex - 1][yIndex].x)
        / (dx * dx);
    const auto d2vxdy = (v[xIndex][yIndex + 1].x - (2 * v[xIndex][yIndex].x) + v[xIndex][yIndex - 1].x)
        / (dy * dy);

    const auto d2vydx = (v[xIndex + 1][yIndex].y - (2 * v[xIndex][yIndex].y) + v[xIndex - 1][yIndex].y)
        / (dx * dx);
    const auto d2vydy = (v[xIndex][yIndex + 1].y - (2 * v[xIndex][yIndex].y) + v[xIndex][yIndex - 1].y)
        / (dy * dy);

    return Vector2d(d2vxdx + d2vxdy, d2vydx + d2vydy);
}