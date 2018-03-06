#pragma once

#include <vector>
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "Render.h"

class Step6NonlinearConvection2D
{
public:
    double fixedTimeStep;

    Step6NonlinearConvection2D()
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

        const auto sigma = 0.2;
        fixedTimeStep = sigma * dx;

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
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        // update U
        auto newUHeights = uHeights;
        for(size_t i = 1; i < uHeights.size(); i++)
        {
            for(size_t j = 1; j < uHeights[i].size(); j++)
            {
                newUHeights[i][j] = uHeights[i][j]
                    - (((uHeights[i][j] * scaledDt) / dx) * (uHeights[i][j] - uHeights[i - 1][j]))
                    - (((vHeights[i][j] * scaledDt) / dy) * (uHeights[i][j] - uHeights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newUHeights);
        uHeights = newUHeights;

        // update V
        auto newVHeights = vHeights;
        for(size_t i = 1; i < vHeights.size(); i++)
        {
            for(size_t j = 1; j < vHeights[i].size(); j++)
            {
                newVHeights[i][j] = vHeights[i][j]
                    - (((uHeights[i][j] * scaledDt) / dx) * (vHeights[i][j] - vHeights[i - 1][j]))
                    - (((vHeights[i][j] * scaledDt) / dy) * (vHeights[i][j] - vHeights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newVHeights);
        vHeights = newVHeights;
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
    const int numX = 81;
    const int numY = 81;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    double dx;
    double dy;
    std::vector<std::vector<double>> uHeights;
    std::vector<std::vector<double>> vHeights;
    SDL_Texture* uHeightMap;
    SDL_Texture* vHeightMap;
};