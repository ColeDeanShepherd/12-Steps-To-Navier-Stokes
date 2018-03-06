#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"

class Step4BurgersEquation1D
{
public:
    double fixedTimeStep;

    Step4BurgersEquation1D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2 * M_PI;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 10;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        nu = 0.07;
        fixedTimeStep = dx * nu;

        // apply initial condition
        const auto nu = 0.3; // viscosity
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);

            const auto xMinus2Pi = x - (2 * M_PI);
            const auto e1 = -((x * x) / (4 * nu));
            const auto e2 = -(pow(xMinus2Pi, 2) / (4 * nu));
            const auto phi = exp(e1) + exp(e2);
            const auto dPhiDx = -(1.0 / (2 * nu)) * ((x * exp(e1)) + (xMinus2Pi * exp(e2)));
            ys[i] = -((2 * nu) / phi) * dPhiDx + 4;
        }
    }
    void update(const double dt)
    {
        const auto timeScale = 1.0 / 4;
        const auto scaledDt = timeScale * dt;
        auto newYs = ys;

        for(size_t i = 1; i < ys.size(); i++)
        {
            const auto iPlus1 = (i + 1) % ys.size();

            newYs[i] = ys[i]
                - (((ys[i] * scaledDt) / dx) * (ys[i] - ys[i - 1]))
                + (((nu * scaledDt) / (dx * dx)) * (ys[iPlus1] - (2 * ys[i]) + ys[i - 1]));
        }

        // apply boundary conditions
        newYs[0] = newYs[newYs.size() - 1];

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 101;
    double dx;
    double nu; // viscosity
    std::vector<double> ys;
};