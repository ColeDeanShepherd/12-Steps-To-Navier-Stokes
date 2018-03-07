#include <vector>
#include <SDL.h>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "LinearAlgebra.h"
#include "Render.h"
#include "Step12.h"

Step12ChannelFlow::Step12ChannelFlow()
{
    title = "Step 12: Channel Flow";

    fixedTimeStep = 1.0 / 60.0;

    graphMetrics.width = WINDOW_WIDTH - 20;
    graphMetrics.height = WINDOW_HEIGHT - 20;
    graphMetrics.pos = Vector2d(10, 10);
    graphMetrics.minX = 0;
    graphMetrics.maxX = 2;
    graphMetrics.minY = 0;
    graphMetrics.maxY = 2;

    dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
    dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

    p = zeros(numX + 2, numY);
    v = std::vector<std::vector<Vector2d>>(numX + 2, std::vector<Vector2d>(numY, Vector2d(0, 0)));
    f = std::vector<std::vector<Vector2d>>(numX + 2, std::vector<Vector2d>(numY, Vector2d(1, 0)));

    heightMap = NULL;
}

void Step12ChannelFlow::applyPBoundaryConditions()
{
    // dp/dy = 0 at y = 0, 2
    for(size_t i = 1; i < (numX + 1); i++)
    {
        p[i][0] = p[i][1];
        p[i][numY - 1] = p[i][numY - 2];
    }
}
void Step12ChannelFlow::updateP(const double dt)
{
    for(size_t iter = 0; iter < numPIterations; iter++)
    {
        auto b = getBForIncompressibleNavierStokes(v, rho, numX + 2, numY, dx, dy, dt);
        for(size_t i = 0; i < numY; i++)
        {
            b[0][i] = b[numX][i];
            b[numX + 1][i] = b[1][i];
        }

        p = iteratePoissonsEquation(p, b, numX + 2, numY, dx, dy);
        applyPBoundaryConditions();
    }
}

void Step12ChannelFlow::applyFlowVelocityBoundaryConditions()
{
    // v = 0 at y = 0, 2
    for(size_t i = 1; i < (numX + 1); i++)
    {
        v[i][0] = Vector2d(0, 0); // v = 0 at y = 0
        v[i][numY - 1] = Vector2d(0, 0); // v = 0 at y = 2
    }
}
void Step12ChannelFlow::updateFlowVelocity(const double dt)
{
    auto newFv = v;
    for(size_t i = 1; i < (numX + 1); i++)
    {
        for(size_t j = 1; j < (numY - 1); j++)
        {
            const auto divergenceOfV = divergence1stOrderBackwardDiff(v, i, j, dx, dy);
            const auto gradientOfP = gradient1stOrderCentralDiff(p, i, j, dx, dy);
            const auto laplacianOfV = laplacian2ndOrderCentralDiff(v, i, j, dx, dy);

            const auto dfvdt = -(divergenceOfV * v[i][j])
                - ((1.0 / rho) * gradientOfP)
                + (nu * laplacianOfV)
                + f[i][j];

            newFv[i][j] = v[i][j] + (dt * (dfvdt));
        }
    }

    v = newFv;

    applyFlowVelocityBoundaryConditions();
}

void Step12ChannelFlow::update(const double dt)
{
    const auto scaledDt = timeScale * dt;

    for(size_t i = 0; i < numY; i++)
    {
        p[0][i] = p[numX][i];
        p[numX + 1][i] = p[1][i];
    }

    for(size_t i = 0; i < numY; i++)
    {
        v[0][i] = v[numX][i];
        v[numX + 1][i] = v[1][i];
    }

    updateP(scaledDt);
    updateFlowVelocity(scaledDt);
}

void Step12ChannelFlow::draw(SDL_Renderer* renderer)
{
    if(heightMap == NULL)
    {
        heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    auto visibleP = std::vector<std::vector<double>>(numX, std::vector<double>(numY));
    for(size_t i = 0; i < numX; i++)
    {
        for(size_t j = 0; j < numY; j++)
        {
            visibleP[i][j] = p[1 + i][j];
        }
    }

    updateHeightmap(heightMap, numX, numY, visibleP, minP, maxP);

    SDL_Rect graphRect;
    graphRect.w = (int)graphMetrics.width;
    graphRect.h = (int)graphMetrics.height;
    graphRect.x = (int)graphMetrics.pos.x;
    graphRect.y = (int)graphMetrics.pos.y;

    SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);

    auto visibleV = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY));
    for(size_t i = 0; i < numX; i++)
    {
        for(size_t j = 0; j < numY; j++)
        {
            visibleV[i][j] = v[1 + i][j];
        }
    }

    const auto pixelHeight = (graphMetrics.maxY - graphMetrics.minY) / numY;
    const auto maxVNorm = (3.0 / 4.0) * pixelHeight;
    renderVectorField(renderer, visibleV, graphMetrics, numX, numY, 0.01, maxVNorm);
}