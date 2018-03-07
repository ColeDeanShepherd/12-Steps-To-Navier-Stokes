#include <vector>
#include <SDL.h>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step9.h"

Step9LaplaceEquation2D::Step9LaplaceEquation2D()
{
    title = "Step 9: Laplace Equation";

    fixedTimeStep = 1.0 / 60;

    graphMetrics.width = WINDOW_WIDTH - 20;
    graphMetrics.height = WINDOW_HEIGHT - 20;
    graphMetrics.pos = Vector2d(10, 10);
    graphMetrics.minX = 0;
    graphMetrics.maxX = 2;
    graphMetrics.minY = 0;
    graphMetrics.maxY = 1;

    dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
    dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

    applyInitialConditions();
    heightMap = NULL;
}
void Step9LaplaceEquation2D::applyInitialConditions()
{
    p = std::vector<std::vector<double>>(numX, std::vector<double>(numY));

    // apply initial condition
    for(size_t i = 0; i < numX; i++)
    {
        for(size_t j = 0; j < numY; j++)
        {
            p[i][j] = 0;
        }
    }
}
void Step9LaplaceEquation2D::applyBoundaryConditions(std::vector<std::vector<double>>& p)
{
    const double bCVal = 1;

    for(size_t i = 0; i < numY; i++)
    {
        p[0][i] = 0; // p = 0 at x = 0

        const auto y = graphMetrics.minY + (i * dy);
        p[numX - 1][i] = y; // p = y at x = 2
    }

    for(size_t i = 0; i < numX; i++)
    {
        p[i][0] = p[i][1]; // dp/dy = 0 at y = 0
        p[i][numY - 1] = p[i][numY - 2]; // dp/dy = 0 at y = 1
    }
}
void Step9LaplaceEquation2D::update(const double dt)
{
    const auto dx2 = dx * dx;
    const auto dy2 = dy * dy;

    auto newP = p;
    for(size_t i = 1; i < (p.size() - 1); i++)
    {
        for(size_t j = 1; j < (p[i].size() - 1); j++)
        {
            const auto numerator = (dy2 * (p[i + 1][j] + p[i - 1][j])) + (dx2 * (p[i][j + 1] + p[i][j - 1]));
            const auto denominator = 2 * (dx2 + dy2);
            newP[i][j] = numerator / denominator;
        }
    }

    applyBoundaryConditions(newP);
    p = newP;
}
void Step9LaplaceEquation2D::draw(SDL_Renderer* renderer)
{
    if(heightMap == NULL)
    {
        heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    updateHeightmap(heightMap, numX, numY, p, 0, graphMaxHeight);

    SDL_Rect uGraphRect;
    uGraphRect.w = (int)graphMetrics.width;
    uGraphRect.h = (int)graphMetrics.height;
    uGraphRect.x = (int)graphMetrics.pos.x;
    uGraphRect.y = (int)graphMetrics.pos.y;

    SDL_RenderCopy(renderer, heightMap, NULL, &uGraphRect);
}