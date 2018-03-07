#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step.h"

class Step3Diffusion1D : public Step
{
public:
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

        fixedTimeStep = sigma * (dx * dx);

        u = std::vector<double>(numPoints);

        // apply initial condition
        for(size_t i = 0; i < numPoints; i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            u[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;
        auto newU = u;

        for(size_t i = 1; i < numPoints - 1; i++)
        {
            const auto d2udx2 = gradient2ndOrderCentralDiff(u, i, dx);
            const auto dudt = nu * d2udx2;
            newU[i] = u[i] + (scaledDt * dudt);
        }

        u = newU;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, u);
    }
private:
    const int numPoints = 41;
    const double sigma = 0.2;
    const double nu = 0.3; // viscosity
    const double timeScale = 1.0 / 10;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};