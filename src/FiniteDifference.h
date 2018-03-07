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
double laplacian2ndOrderCentralDiff(
    const std::vector<std::vector<double>>& f,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto d2fdx2 = (f[xIndex + 1][yIndex] - (2 * f[xIndex][yIndex]) + f[xIndex - 1][yIndex])
        / (dx * dx);
    const auto d2fdy2 = (f[xIndex][yIndex + 1] - (2 * f[xIndex][yIndex]) + f[xIndex][yIndex - 1])
        / (dy * dy);

    return d2fdx2 + d2fdy2;
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

std::vector<std::vector<double>> iteratePoissonsEquation(
    const std::vector<std::vector<double>>& p, const std::vector<std::vector<double>>& b,
    const size_t numX, const size_t numY, const double dx, const double dy)
{
    const auto dx2 = dx * dx;
    const auto dy2 = dy * dy;

    auto newP = p;
    for(size_t i = 1; i < (numX - 1); i++)
    {
        for(size_t j = 1; j < (numY - 1); j++)
        {
            const auto term1Numerator =
                (dy2 * (p[i + 1][j] + p[i - 1][j])) +
                (dx2 * (p[i][j + 1] + p[i][j - 1]));
            const auto termDenominator = 2 * (dx2 + dy2);
            const auto term1 = term1Numerator / termDenominator;
            const auto term2 = ((dx2 * dy2) / termDenominator) * b[i][j];

            newP[i][j] = term1 - term2;
        }
    }

    return newP;
}