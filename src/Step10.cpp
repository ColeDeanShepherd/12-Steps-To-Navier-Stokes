#include <vector>
#include <SDL.h>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "LinearAlgebra.h"
#include "Render.h"
#include "Step10.h"

Step10PoissonEquation2D::Step10PoissonEquation2D()
{
    title = "Step 10: Poisson Equation";

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

    // init p
    p = zeros(numX, numY);

    // init b
    b = zeros(numX, numY);
    b[(int)((double)numX / 4)][(int)((double)numY / 4)] = 100;
    b[(int)(3 * (double)numX / 4)][(int)(3 * (double)numY / 4)] = -100;

    heightMap = NULL;
}
void Step10PoissonEquation2D::applyBoundaryConditions()
{
    const double bCVal = 1;

    for(size_t i = 0; i < numY; i++)
    {
        p[0][i] = 0; // p = 0 at x = 0
        p[numX - 1][i] = 0; // p = 0 at x = 2
    }

    for(size_t i = 0; i < numX; i++)
    {
        p[i][0] = 0; // p = 0 at y = 0
        p[i][numY - 1] = 0; // p = 0 at y = 1
    }
}
void Step10PoissonEquation2D::update(const double dt)
{
    p = iteratePoissonsEquation(p, b, numX, numY, dx, dy);
    applyBoundaryConditions();
}
void Step10PoissonEquation2D::draw(SDL_Renderer* renderer)
{
    if(heightMap == NULL)
    {
        heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    updateHeightmap(heightMap, numX, numY, p, graphMinHeight, graphMaxHeight);

    SDL_Rect uGraphRect;
    uGraphRect.w = (int)graphMetrics.width;
    uGraphRect.h = (int)graphMetrics.height;
    uGraphRect.x = (int)graphMetrics.pos.x;
    uGraphRect.y = (int)graphMetrics.pos.y;

    SDL_RenderCopy(renderer, heightMap, NULL, &uGraphRect);
}