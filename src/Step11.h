#pragma once

#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "LinearAlgebra.h"
#include "Render.h"
#include "Step.h"

class Step11CavityFlow : public Step
{
public:
    Step11CavityFlow()
    {
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

        b = zeros(numX, numY);
        p = zeros(numX, numY);
        fv = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY, Vector2d(0, 0)));

        heightMap = NULL;
    }

    void buildUpB(const double dt)
    {
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                const auto term1 = (1.0 / dt) * divergence1stOrderBackwardDiff(fv, i, j, dx, dy);
                const auto term2 = pow((fv[i + 1][j].x - fv[i - 1][j].x) / (2 * dx), 2);
                const auto term3 = 2 *
                    ((fv[i][j + 1].x - fv[i][j - 1].x) / (2 * dy)) *
                    ((fv[i + 1][j].y - fv[i - 1][j].y) / (2 * dx));
                const auto term4 = pow((fv[i][j + 1].y - fv[i][j - 1].y) / (2 * dy), 2);

                b[i][j] = rho * (term1 - term2 - term3 - term4);
            }
        }
    }
    void applyPBoundaryConditions()
    {
        // dp/dy = 0 at y = 0
        for(size_t i = 0; i < numX; i++)
        {
            p[i][0] = p[i][1];
        }

        // p = 0 at y = 2
        for(size_t i = 0; i < numX; i++)
        {
            p[i][numY - 1] = 0;
        }

        // dp/dx = 0 at x = 0, 2
        for(size_t i = 0; i < numX; i++)
        {
            // dp/dx = 0 at x = 0
            p[0][i] = p[1][i];

            // dp/dx = 0 at x = 2
            p[numX - 1][i] = p[numX - 2][i];
        }
    }
    void updateP(const double dt)
    {
        for(size_t iter = 0; iter < numPIterations; iter++)
        {
            buildUpB(dt);
            p = iteratePoissonsEquation(p, b, numX, numY, dx, dy);
            applyPBoundaryConditions();
        }
    }

    void applyFlowVelocityBoundaryConditions()
    {
        // v = 0 at boundaries
        for(size_t i = 0; i < numX; i++)
        {
            fv[0][i] = Vector2d(0, 0); // v = 0 at x = 0
            fv[numX - 1][i] = Vector2d(0, 0); // v = 0 at x = 2
            fv[i][0] = Vector2d(0, 0); // v = 0 at y = 0
            fv[i][numY - 1] = Vector2d(0, 0); // v = 0 at y = 2
        }

        // vx = 1 at y = 2
        for(size_t i = 0; i < numX; i++)
        {
            fv[i][numY - 1].x = 1;
        }
    }
    void updateFlowVelocity(const double dt)
    {
        auto newFv = fv;
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                const auto divergenceOfV = divergence1stOrderBackwardDiff(fv, i, j, dx, dy);
                const auto gradientOfP = gradient1stOrderCentralDiff(p, i, j, dx, dy);
                const auto laplacianOfV = laplacian2ndOrderCentralDiff(fv, i, j, dx, dy);

                const auto dfvdt = -(divergenceOfV * fv[i][j])
                    - ((1.0 / rho) * gradientOfP)
                    + (nu * laplacianOfV);

                newFv[i][j] = fv[i][j] + (dt * (dfvdt));
            }
        }

        fv = newFv;

        applyFlowVelocityBoundaryConditions();
    }

    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        updateP(scaledDt);
        updateFlowVelocity(scaledDt);
    }

    void draw(SDL_Renderer* renderer)
    {
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap(heightMap, numX, numY, p, minP, maxP);

        SDL_Rect graphRect;
        graphRect.w = (int)graphMetrics.width;
        graphRect.h = (int)graphMetrics.height;
        graphRect.x = (int)graphMetrics.pos.x;
        graphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);

        renderVectorField(renderer, fv, graphMetrics, numX, numY, 1);
    }
private:
    const int numX = 41;
    const int numY = 41;
    const double timeScale = 0.006;
    const int numPIterations = 50;
    const double c = 1;
    const double rho = 1;
    const double nu = 0.1;
    const double minP = -4;
    const double maxP = 4;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> fv;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    SDL_Texture* heightMap;
};