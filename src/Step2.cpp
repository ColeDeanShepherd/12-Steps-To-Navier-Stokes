#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step2.h"

Step2NonlinearConvection1D::Step2NonlinearConvection1D()
{
    title = "Step 2: Nonlinear Convection";

    fixedTimeStep = 1.0 / 60.0;

    graphMetrics.width = WINDOW_WIDTH - 20;
    graphMetrics.height = WINDOW_HEIGHT - 20;
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
void Step2NonlinearConvection1D::update(const double dt)
{
    const auto timeScale = 1.0 / 4;
    const auto scaledDt = timeScale * dt;
    auto newU = u;

    for(size_t i = 1; i < numPoints; i++)
    {
        const auto dudx = gradient1stOrderBackwardDiff(u, i, dx);
        const auto dudt = -(u[i] * dudx);
        newU[i] = u[i] + scaledDt * dudt;
    }

    u = newU;
}
void Step2NonlinearConvection1D::draw(SDL_Renderer* renderer)
{
    renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, u);
}