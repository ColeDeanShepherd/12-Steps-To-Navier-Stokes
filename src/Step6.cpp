#include <vector>
#include <SDL.h>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step6.h"

Step6NonlinearConvection2D::Step6NonlinearConvection2D(const int windowWidth, const int windowHeight)
{
    title = "Step 6: 2D Nonlinear Convection";

    graphMetrics.width = windowWidth - 20;
    graphMetrics.height = windowHeight - 20;
    graphMetrics.pos = Vector2d(10, 10);
    graphMetrics.minX = 0;
    graphMetrics.maxX = 2;
    graphMetrics.minY = 0;
    graphMetrics.maxY = 2;

    dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
    dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

    fixedTimeStep = sigma * dx;

    applyInitialCondition();

    vxHeightMap = NULL;
    vyHeightMap = NULL;
}
void Step6NonlinearConvection2D::applyInitialCondition()
{
    v = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY));

    for(size_t i = 0; i < numX; i++)
    {
        const auto x = graphMetrics.minX + (i * dx);

        for(size_t j = 0; j < numY; j++)
        {
            const auto y = graphMetrics.minY + (j * dy);
            const auto value = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;

            v[i][j] = Vector2d(value, value);
        }
    }
}
void Step6NonlinearConvection2D::applyBoundaryConditions()
{
    const double bCVal = 1;

    for(size_t i = 0; i < numY; i++)
    {
        v[0][i] = Vector2d(bCVal, bCVal); // left
        v[numX - 1][i] = Vector2d(bCVal, bCVal); // right
    }

    for(size_t i = 0; i < numX; i++)
    {
        v[i][0] = Vector2d(bCVal, bCVal); // bottom
        v[i][numY - 1] = Vector2d(bCVal, bCVal); // top
    }
}
void Step6NonlinearConvection2D::update(const double dt)
{
    const auto scaledDt = timeScale * dt;

    // update V
    auto newV = v;
    for(size_t i = 1; i < numX; i++)
    {
        for(size_t j = 1; j < numY; j++)
        {
            const auto gradVx = Vector2d(
                (v[i][j].x - v[i - 1][j].x) / dx,
                (v[i][j].x - v[i][j - 1].x) / dy
            );
            const auto gradVy = Vector2d(
                (v[i][j].y - v[i - 1][j].y) / dx,
                (v[i][j].y - v[i][j - 1].y) / dy
            );

            const auto dvdt = Vector2d(
                -dot(v[i][j], gradVx),
                -dot(v[i][j], gradVy)
            );

            newV[i][j] = v[i][j] + (scaledDt * dvdt);
        }
    }

    v = newV;
    applyBoundaryConditions();
}
void Step6NonlinearConvection2D::draw(SDL_Renderer* renderer)
{
    if(vxHeightMap == NULL)
    {
        vxHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    if(vyHeightMap == NULL)
    {
        vyHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
    }

    std::vector<std::vector<double>> vx(numX, std::vector<double>(numY, 0));
    std::vector<std::vector<double>> vy(numX, std::vector<double>(numY, 0));

    for(size_t i = 0; i < numX; i++)
    {
        for(size_t j = 0; j < numY; j++)
        {
            vx[i][j] = v[i][j].x;
            vy[i][j] = v[i][j].y;
        }
    }

    updateHeightmap(vxHeightMap, numX, numY, vx, 0, graphMaxHeight);
    updateHeightmap(vyHeightMap, numX, numY, vy, 0, graphMaxHeight);

    SDL_Rect vxGraphRect;
    vxGraphRect.w = (int)((graphMetrics.width / 2) - 5);
    vxGraphRect.h = (int)graphMetrics.height;
    vxGraphRect.x = (int)graphMetrics.pos.x;
    vxGraphRect.y = (int)graphMetrics.pos.y;

    SDL_RenderCopy(renderer, vxHeightMap, NULL, &vxGraphRect);

    SDL_Rect vyGraphRect;
    vyGraphRect.w = vxGraphRect.w;
    vyGraphRect.h = (int)graphMetrics.height;
    vyGraphRect.x = (int)(vxGraphRect.x + vxGraphRect.w + 10);
    vyGraphRect.y = (int)graphMetrics.pos.y;

    SDL_RenderCopy(renderer, vyHeightMap, NULL, &vyGraphRect);
}