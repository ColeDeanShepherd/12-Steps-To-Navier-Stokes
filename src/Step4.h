#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step.h"

class Step4BurgersEquation1D : public Step
{
public:
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

        fixedTimeStep = dx * nu;

        u.resize(numPoints);

        // apply initial condition
        const auto icNu = 0.3; // viscosity
        for(size_t i = 0; i < numPoints; i++)
        {
            const auto x = graphMetrics.minX + (i * dx);

            const auto xMinus2Pi = x - (2 * M_PI);
            const auto e1 = -((x * x) / (4 * icNu));
            const auto e2 = -(pow(xMinus2Pi, 2) / (4 * icNu));
            const auto phi = exp(e1) + exp(e2);
            const auto dPhiDx = -(1.0 / (2 * icNu)) * ((x * exp(e1)) + (xMinus2Pi * exp(e2)));
            u[i] = -((2 * icNu) / phi) * dPhiDx + 4;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;
        auto newU = u;

        for(size_t i = 1; i < numPoints; i++)
        {
            const auto iPlus1 = (i + 1) % numPoints;
            
            const auto dudx = gradient1stOrderBackwardDiff(u, i, dx);
            const auto d2udx2 = (u[iPlus1] - (2 * u[i]) + u[i - 1]) / (dx * dx);
            const auto dudt = (nu * d2udx2) - (u[i] * dudx);

            newU[i] = u[i] + (scaledDt * dudt);
        }

        // apply boundary conditions
        newU[0] = newU[newU.size() - 1];

        u = newU;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, u);
    }
private:
    const int numPoints = 101;
    const double nu = 0.07; // viscosity
    const double timeScale = 1.0 / 4;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};