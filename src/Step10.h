#pragma once

#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "LinearAlgebra.h"
#include "Render.h"

class Step10PoissonEquation2D
{
public:
    double fixedTimeStep;

    Step10PoissonEquation2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 1;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        fixedTimeStep = 1.0 / 60;

        // init p
        p = zeros(numY, numX);

        // init b
        b = zeros(numY, numX);
        b[(int)((double)numY / 4)][(int)((double)numX / 4)] = 100;
        b[(int)(3 * (double)numY / 4)][(int)(3 * (double)numX / 4)] = -100;

        heightMap = NULL;
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& p)
    {
        const double bCVal = 1;

        // p = 0 at x = 0
        for(size_t i = 0; i < numY; i++)
        {
            p[i][0] = 0;
        }

        // p = 0 at x = 2
        for(size_t i = 0; i < numY; i++)
        {
            p[i][numX - 1] = 0;
        }

        // p = 0 at y = 0
        for(size_t i = 0; i < numY; i++)
        {
            p[0][i] = 0;
        }

        // p = 0 at y = 1
        for(size_t i = 0; i < numY; i++)
        {
            p[numY - 1][i] = 0;
        }
    }
    void update(std::vector<std::vector<double>>& p, const double dt)
    {
        const auto dx2 = dx * dx;
        const auto dy2 = dy * dy;

        auto newP = p;
        for(size_t i = 1; i < (p.size() - 1); i++)
        {
            for(size_t j = 1; j < (p[i].size() - 1); j++)
            {
                const auto numerator =
                    (dy2 * (p[i + 1][j] + p[i - 1][j])) +
                    (dx2 * (p[i][j + 1] + p[i][j - 1])) -
                    (b[i][j] * dx2 * dy2);
                const auto denominator = 2 * (dx2 + dy2);
                newP[i][j] = numerator / denominator;
            }
        }

        applyBoundaryConditions(newP);
        p = newP;
    }
    void update(const double dt)
    {
        update(p, dt);
    }
    void updateHeightmap(const std::vector<std::vector<double>>& heights, SDL_Texture* heightMap)
    {
        unsigned char* pixelBytes;
        int pitch;
        const auto failedLocking = SDL_LockTexture(heightMap, NULL, (void**)&pixelBytes, &pitch) != 0;

        if(failedLocking)
        {
            throw new std::exception("Failed locking the texture.");
        }

        for(auto rowIndex = 0; rowIndex < numY; rowIndex++)
        {
            for(auto columnIndex = 0; columnIndex < numX; columnIndex++)
            {
                const auto pixelBytesOffset = (rowIndex * pitch) + (4 * columnIndex);
                const auto pixelValuePercent =
                    (heights[rowIndex][columnIndex] - graphMinHeight) /
                    (graphMaxHeight - graphMinHeight);
                const auto pixelValue = (unsigned char)(255 * pixelValuePercent);

                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void draw(SDL_Renderer* renderer)
    {
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap(p, heightMap);

        SDL_Rect uGraphRect;
        uGraphRect.w = (int)graphMetrics.width;
        uGraphRect.h = (int)graphMetrics.height;
        uGraphRect.x = (int)graphMetrics.pos.x;
        uGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, heightMap, NULL, &uGraphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 50;
    const int numY = 50;
    const double graphMinHeight = -0.06;
    const double graphMaxHeight = 0.06;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    SDL_Texture* heightMap;
};