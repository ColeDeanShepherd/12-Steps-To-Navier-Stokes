#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"

class Step5LinearConvection2D
{
public:
    double fixedTimeStep;

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

        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        const auto sigma = 0.2;
        fixedTimeStep = sigma * dx;

        applyInitialConditions();
        heightMap = NULL;
    }
    void applyInitialConditions()
    {
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            const auto y = graphMetrics.minY + (rowIndex * dy);

            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                const auto x = graphMetrics.minX + (columnIndex * dx);
                heights[rowIndex][columnIndex] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
            }
        }
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& newHeights)
    {
        const double bCVal = 1;

        // left
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][0] = bCVal;
        }
        // right
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][numX - 1] = bCVal;
        }
        // top
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[0][i] = bCVal;
        }
        // bottom
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[numY - 1][i] = bCVal;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        auto newHeights = heights;

        // update
        for(size_t i = 1; i < heights.size(); i++)
        {
            for(size_t j = 1; j < heights[i].size(); j++)
            {
                newHeights[i][j] = heights[i][j]
                    - (((c * scaledDt) / dx) * (heights[i][j] - heights[i - 1][j]))
                    - (((c * scaledDt) / dy) * (heights[i][j] - heights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newHeights);

        heights = newHeights;
    }
    void updateHeightmap()
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
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap();

        SDL_Rect graphRect;
        graphRect.w = (int)graphMetrics.width;
        graphRect.h = (int)graphMetrics.height;
        graphRect.x = (int)graphMetrics.pos.x;
        graphRect.y = (int)graphMetrics.pos.y;
        SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 81;
    const int numY = 81;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    double dx;
    double dy;
    std::vector<std::vector<double>> heights;
    SDL_Texture* heightMap;
};