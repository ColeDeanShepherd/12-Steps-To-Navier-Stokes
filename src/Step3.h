#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"

class Step3Diffusion1D
{
public:
    double fixedTimeStep;

    Step3Diffusion1D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        const auto sigma = 0.2f;
        fixedTimeStep = sigma * (dx * dx);

        // apply initial condition
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            ys[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto timeScale = 1.0 / 10;
        const auto nu = 0.3f; // viscosity
        auto newYs = ys;

        for(size_t i = 1; i < ys.size() - 1; i++)
        {
            newYs[i] = ys[i] + (nu * ((timeScale * dt) / (dx * dx)) * (ys[i + 1] - (2 * ys[i]) + ys[i - 1]));
        }

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 41;
    double dx;
    std::vector<double> ys;
};