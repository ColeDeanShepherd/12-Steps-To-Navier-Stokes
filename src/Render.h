#pragma once

#include <vector>
#include "Core.h"
#include "GraphMetrics.h"

struct SDL_Renderer;
struct SDL_Texture;

void renderLineGraph(SDL_Renderer* renderer, const GraphMetrics& graphMetrics, const double x0, const double dx, const std::vector<double>& ys);

void updateHeightmap(
    SDL_Texture* heightMap, const size_t width, const size_t height,
    const std::vector<std::vector<double>>& values, const double minValue, const double maxValue);

void renderVectorField(
    SDL_Renderer* renderer, const std::vector<std::vector<Vector2d>>& v, const GraphMetrics& graphMetrics,
    const size_t numX, const size_t numY, const double vScale, const double maxVNorm = -1);