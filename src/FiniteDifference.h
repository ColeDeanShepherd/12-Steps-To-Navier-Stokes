#pragma once

#include <vector>
#include "Vector2d.h"

double gradient1stOrderBackwardDiff(const std::vector<double>& f, const size_t i, const double dx);
double gradient2ndOrderCentralDiff(const std::vector<double>& f, const size_t i, const double dx);
Vector2d gradient1stOrderCentralDiff(
    const std::vector<std::vector<double>>& f,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
);
double divergence1stOrderBackwardDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
);
double laplacian2ndOrderCentralDiff(
    const std::vector<std::vector<double>>& f,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
);
Vector2d laplacian2ndOrderCentralDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
);

std::vector<std::vector<double>> iteratePoissonsEquation(
    const std::vector<std::vector<double>>& p, const std::vector<std::vector<double>>& b,
    const size_t numX, const size_t numY, const double dx, const double dy);

std::vector<std::vector<double>> getBForIncompressibleNavierStokes(
    const std::vector<std::vector<Vector2d>>& v, const double rho,
    const size_t numX, const size_t numY, const double dx, const double dy, const double dt);