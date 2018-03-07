#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"
#include "Step.h"

class Step5LinearConvection2D : public Step
{
public:
    Step5LinearConvection2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        fixedTimeStep = sigma * dx;

        u = std::vector<std::vector<double>>(numX, std::vector<double>(numY));

        applyInitialConditions();
        heightMap = NULL;
    }
    void applyInitialConditions()
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
    void applyBoundaryConditions(std::vector<std::vector<double>>& newHeights)
    {
        const double bCVal = 1;

        for(size_t i = 0; i < numY; i++)
        {
            newHeights[0][i] = bCVal; // left
            newHeights[numX - 1][i] = bCVal; // right
        }
        
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[i][0] = bCVal; // bottom
            newHeights[i][numY - 1] = bCVal; // top
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        auto newU = u;

        for(size_t i = 1; i < numX; i++)
        {
            for(size_t j = 1; j < numY; j++)
            {
                const auto dudx = (u[i][j] - u[i - 1][j]) / dx;
                const auto dudy = (u[i][j] - u[i][j - 1]) / dy;
                const auto dudt = -(c * dudx) - (c * dudy);

                newU[i][j] = u[i][j] + (scaledDt * dudt);
            }
        }

        applyBoundaryConditions(newU);

        u = newU;
    }
    void draw(SDL_Renderer* renderer)
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
private:
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double sigma = 0.2;
    const int numX = 81;
    const int numY = 81;
    
    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> u;
    SDL_Texture* heightMap;
};