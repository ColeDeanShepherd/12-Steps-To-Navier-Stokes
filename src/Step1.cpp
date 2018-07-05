#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step1.h"

Step1LinearConvection1D::Step1LinearConvection1D(const int windowWidth, const int windowHeight)
{
    title = "Step 1: Linear Convection";

    fixedTimeStep = 1.0 / 60.0;

    graphMetrics.width = windowWidth - 20;
    graphMetrics.height = windowHeight - 20;
    graphMetrics.pos = Vector2d(10, 10);
    graphMetrics.minX = 0;
    graphMetrics.maxX = 2;
    graphMetrics.minY = 0;
    graphMetrics.maxY = 2;

    dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

    u = std::vector<double>(numPoints);

    // apply initial condition
    for(size_t i = 0; i < numPoints; i++)
    {
        const auto x = graphMetrics.minX + (i * dx);
        u[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
    }
}
void Step1LinearConvection1D::update(const double dt)
{
    auto newU = u;

    for(size_t i = 1; i < numPoints; i++)
    {
        const auto dudx = gradient1stOrderBackwardDiff(u, i, dx);
        const auto dudt = -(c * dudx);
        newU[i] = u[i] + (dt * dudt);
    }

    u = newU;
}
void Step1LinearConvection1D::draw(SDL_Renderer* renderer)
{
    renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, u);
}