#pragma once

#include <vector>
#include "Core.h"
#include "Vector2d.h"
#include "GraphMetrics.h"
#include "FiniteDifference.h"
#include "LinearAlgebra.h"
#include "Render.h"

class Step12ChannelFlow
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step12ChannelFlow()
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

        b = zeros(numX, numY);
        p = zeros(numX, numY);
        fv = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY, Vector2d(0, 0)));
        f = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY, Vector2d(1, 0)));

        heightMap = NULL;
    }

    void buildUpB(const double dt)
    {
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                b[i][j] = getB(i, j, dt);
            }
        }
    }
    double getB(const size_t i, const size_t j, const double dt)
    {
        const auto term1 = (1.0 / dt) * divergence1stOrderBackwardDiff(fv, i, j, dx, dy);
        const auto term2 = pow((fv[i + 1][j].x - fv[i - 1][j].x) / (2 * dx), 2);
        const auto term3 = 2 *
            ((fv[i][j + 1].x - fv[i][j - 1].x) / (2 * dy)) *
            ((fv[i + 1][j].y - fv[i - 1][j].y) / (2 * dx));
        const auto term4 = pow((fv[i][j + 1].y - fv[i][j - 1].y) / (2 * dy), 2);

        return term1 - term2 - term3 - term4;
    }
    void iterateP(const double dt)
    {
        auto newP = p;
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                const auto term1Numerator =
                    ((dy * dy) * (p[i + 1][j] + p[i - 1][j])) +
                    ((dx * dx) * (p[i][j + 1] + p[i][j - 1]));
                const auto termDenominator = 2 * ((dx * dx) + (dy * dy));
                const auto term1 = term1Numerator / termDenominator;
                const auto term2 = ((rho * (dx * dx) * (dy * dy)) / termDenominator) * getB(i, j, dt);

                newP[i][j] = term1 - term2;
            }
        }

        p = newP;

        // dp/dy = 0 at y = 0, 2
        for(size_t i = 0; i < numX; i++)
        {
            p[i][0] = p[i][1];
            p[i][numY - 1] = p[i][numY - 2];
        }
    }
    void updateP(const double dt)
    {
        for(size_t iter = 0; iter < numPIterations; iter++)
        {
            iterateP(dt);
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
                    + (nu * laplacianOfV)
                    + f[i][j];

                newFv[i][j] = fv[i][j] + (dt * (dfvdt));
            }
        }

        fv = newFv;

        // v = 0 at y = 0, 2
        for(size_t i = 0; i < numX; i++)
        {
            fv[i][0] = Vector2d(0, 0); // v = 0 at y = 0
            fv[i][numY - 1] = Vector2d(0, 0); // v = 0 at y = 2
        }
    }

    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        updateP(scaledDt);
        updateFlowVelocity(scaledDt);
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

        // pixels go from top to bottom, left to right
        for(auto pixelRowIndex = 0; pixelRowIndex < numY; pixelRowIndex++)
        {
            for(auto pixelColumnIndex = 0; pixelColumnIndex < numX; pixelColumnIndex++)
            {
                const auto pixelBytesOffset = (pitch * pixelRowIndex) + (4 * pixelColumnIndex);

                const auto pPercent = std::clamp<double>(
                    (p[pixelColumnIndex][(numY - 1) - pixelRowIndex] - minP) / (maxP - minP),
                    0, 1);

                const auto pixelValue = (unsigned char)(255 * pPercent);

                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void drawVelocityVectors(SDL_Renderer* renderer)
    {
        for(size_t i = 0; i < numX; i++)
        {
            for(size_t j = 0; j < numY; j++)
            {
                const auto pixelWidth = (graphMetrics.maxX - graphMetrics.minX) / numX;
                const auto pixelHeight = (graphMetrics.maxY - graphMetrics.minY) / numY;

                const auto arrow = withMaxNorm(fv[i][j], (3.0 / 4.0) * pixelHeight);
                const auto velocityTailGraphPos = Vector2d(
                    graphMetrics.minX + (i * pixelWidth) + (pixelWidth / 2),
                    graphMetrics.minY + (j * pixelHeight) + (pixelHeight / 2));
                const auto velocityHeadGraphPos = add(velocityTailGraphPos, arrow);
                const auto velocityTailPxPos = graphPointToPxPoint(graphMetrics, velocityTailGraphPos);
                const auto velocityHeadPxPos = graphPointToPxPoint(graphMetrics, velocityHeadGraphPos);

                SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                SDL_RenderDrawLine(
                    renderer, velocityTailPxPos.x, velocityTailPxPos.y,
                    velocityHeadPxPos.x, velocityHeadPxPos.y
                );

                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
                SDL_RenderDrawPoint(renderer, velocityHeadPxPos.x, velocityHeadPxPos.y);
            }
        }
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

        // draw velocity vectors
        drawVelocityVectors(renderer);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 41;
    const int numY = 41;
    const double timeScale = 0.006;
    const int numPIterations = 50;
    const double c = 1;
    const double rho = 1;
    const double nu = 0.1;
    const double minP = -4;
    const double maxP = 4;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> fv;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    std::vector<std::vector<Vector2d>> f;
    SDL_Texture* heightMap;
};