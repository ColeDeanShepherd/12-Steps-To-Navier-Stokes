#include <vector>
#include <SDL.h>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step7.h"

Step7Diffusion2D::Step7Diffusion2D()
{
    title = "Step 7: 2D Diffusion";

    graphMetrics.width = WINDOW_WIDTH - 20;
    graphMetrics.height = WINDOW_HEIGHT - 20;
    graphMetrics.pos = Vector2d(10, 10);
    graphMetrics.minX = 0;
    graphMetrics.maxX = 2;
    graphMetrics.minY = 0;
    graphMetrics.maxY = 2;

    dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
    dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

    fixedTimeStep = (sigma * dx * dy) / nu;

    u = std::vector<std::vector<double>>(numX, std::vector<double>(numY));

    applyInitialConditions();
    heightMap = NULL;
}
void Step7Diffusion2D::applyInitialConditions()
{
    for(size_t i = 0; i < numX; i++)
    {
        const auto x = graphMetrics.minX + (i * dx);

        for(size_t j = 0; j < numY; j++)
        {
            const auto y = graphMetrics.minY + (j * dy);
            u[i][j] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
        }
    }
}
void Step7Diffusion2D::applyBoundaryConditions(std::vector<std::vector<double>>& newU)
{
    const double bCVal = 1;

    for(size_t i = 0; i < numY; i++)
    {
        newU[0][i] = bCVal; // left
        newU[numX - 1][i] = bCVal; // right
    }

    for(size_t i = 0; i < numX; i++)
    {
        newU[i][0] = bCVal; // bottom
        newU[i][numY - 1] = bCVal; // top
    }
}
void Step7Diffusion2D::update(const double dt)
{
    const auto scaledDt = timeScale * dt;

    auto newU = u;

    // update
    for(size_t i = 1; i < (numX - 1); i++)
    {
        for(size_t j = 1; j < (numY - 1); j++)
        {
            const auto laplacianU = laplacian2ndOrderCentralDiff(u, i, j, dx, dy);
            const auto dudt = nu * laplacianU;

            newU[i][j] = u[i][j] + (scaledDt * dudt);
        }
    }

    applyBoundaryConditions(newU);

    u = newU;
}
void Step7Diffusion2D::draw(SDL_Renderer* renderer)
{
    if(heightMap == NULL)
    {
        heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    updateHeightmap(heightMap, numX, numY, u, 0, graphMaxHeight);

    SDL_Rect graphRect;
    graphRect.w = (int)graphMetrics.width;
    graphRect.h = (int)graphMetrics.height;
    graphRect.x = (int)graphMetrics.pos.x;
    graphRect.y = (int)graphMetrics.pos.y;
    SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);
}