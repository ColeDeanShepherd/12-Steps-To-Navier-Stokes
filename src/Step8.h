#pragma once

#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"

class Step8BurgersEquation2D
{
public:
    double fixedTimeStep;

    Step8BurgersEquation2D()
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

        const auto sigma = 0.0009;
        fixedTimeStep = (sigma * dx * dy) / nu;

        initForFunction(uHeights, &uHeightMap);
        initForFunction(vHeights, &vHeightMap);
    }
    void initForFunction(std::vector<std::vector<double>>& heights, SDL_Texture** heightMap)
    {
        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        // apply initial condition
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            const auto y = graphMetrics.minY + (rowIndex * dy);

            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                const auto x = graphMetrics.minX + (columnIndex * dx);
                heights[rowIndex][columnIndex] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
            }
        }

        heightMap = NULL;
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& heights)
    {
        const double bCVal = 1;

        // left
        for(size_t i = 0; i < numY; i++)
        {
            heights[i][0] = bCVal;
        }
        // right
        for(size_t i = 0; i < numY; i++)
        {
            heights[i][numX - 1] = bCVal;
        }
        // top
        for(size_t i = 0; i < numX; i++)
        {
            heights[0][i] = bCVal;
        }
        // bottom
        for(size_t i = 0; i < numX; i++)
        {
            heights[numY - 1][i] = bCVal;
        }
    }
    void update(std::vector<std::vector<double>>& heights, const double dt)
    {
        const auto scaledDt = timeScale * dt;
        const auto& u = uHeights;
        const auto& v = vHeights;
        auto w = heights.data();

        auto newHeights = heights;
        for(size_t i = 1; i < (heights.size() - 1); i++)
        {
            for(size_t j = 1; j < (heights[i].size() - 1); j++)
            {
                newHeights[i][j] = w[i][j]
                    - (((u[i][j] * scaledDt) / dx) * (w[i][j] - w[i - 1][j]))
                    - (((v[i][j] * scaledDt) / dy) * (w[i][j] - w[i][j - 1]))
                    + (((nu * scaledDt) / (dx * dx)) * (w[i + 1][j] - (2 * w[i][j]) + w[i - 1][j]))
                    + (((nu * scaledDt) / (dy * dy)) * (w[i][j + 1] - (2 * w[i][j]) + w[i][j - 1]));
            }
        }

        applyBoundaryConditions(newHeights);
        heights = newHeights;
    }
    void update(const double dt)
    {
        update(uHeights, dt);
        update(vHeights, dt);
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
                const auto pixelValuePercent = heights[rowIndex][columnIndex] / graphMaxHeight;
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
        if(uHeightMap == NULL)
        {
            uHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        if(vHeightMap == NULL)
        {
            vHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap(uHeights, uHeightMap);
        updateHeightmap(vHeights, vHeightMap);

        SDL_Rect uGraphRect;
        uGraphRect.w = ((int)graphMetrics.width / 2) - 5;
        uGraphRect.h = (int)graphMetrics.height;
        uGraphRect.x = (int)graphMetrics.pos.x;
        uGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, uHeightMap, NULL, &uGraphRect);

        SDL_Rect vGraphRect;
        vGraphRect.w = ((int)graphMetrics.width / 2) - 5;
        vGraphRect.h = (int)graphMetrics.height;
        vGraphRect.x = uGraphRect.x + uGraphRect.w + 10;
        vGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, vHeightMap, NULL, &vGraphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 41;
    const int numY = 41;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double nu = 0.01;
    double dx;
    double dy;
    std::vector<std::vector<double>> uHeights;
    std::vector<std::vector<double>> vHeights;
    SDL_Texture* uHeightMap;
    SDL_Texture* vHeightMap;
};