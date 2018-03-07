#include <SDL.h>
#include "Render.h"

void renderLineGraph(SDL_Renderer* renderer, const GraphMetrics& graphMetrics, const double x0, const double dx, const std::vector<double>& ys)
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    for(auto i = 0; i < ys.size() - 1; i++)
    {
        const auto graphPoint0 = Vector2d(x0 + (i * dx), ys[i]);
        const auto pixelPoint0 = graphPointToPxPoint(graphMetrics, graphPoint0);

        const auto graphPoint1 = Vector2d(x0 + ((i + 1) * dx), ys[i + 1]);
        const auto pixelPoint1 = graphPointToPxPoint(graphMetrics, graphPoint1);

        SDL_RenderDrawLine(renderer, (int)pixelPoint0.x, (int)pixelPoint0.y, (int)pixelPoint1.x, (int)pixelPoint1.y);
    }
}

void updateHeightmap(
    SDL_Texture* heightMap, const size_t width, const size_t height,
    const std::vector<std::vector<double>>& values, const double minValue, const double maxValue)
{
    unsigned char* pixelBytes;
    int pitch;
    const auto failedLocking = SDL_LockTexture(heightMap, NULL, (void**)&pixelBytes, &pitch) != 0;

    if(failedLocking)
    {
        throw new std::exception("Failed locking the texture.");
    }

    // pixels go from top to bottom, left to right
    for(auto pixelRowIndex = 0; pixelRowIndex < height; pixelRowIndex++)
    {
        for(auto pixelColumnIndex = 0; pixelColumnIndex < width; pixelColumnIndex++)
        {
            const auto pixelBytesOffset = (pitch * pixelRowIndex) + (4 * pixelColumnIndex);

            const auto valuePercent = std::clamp<double>(
                (values[pixelColumnIndex][(height - 1) - pixelRowIndex] - minValue) / (maxValue - minValue),
                0, 1);
            const auto pixelValue = (unsigned char)(255 * valuePercent);

            pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
            pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
            pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
            pixelBytes[pixelBytesOffset + 0] = 255; // A
        }
    }

    SDL_UnlockTexture(heightMap);
}

void renderVectorField(
    SDL_Renderer* renderer, const std::vector<std::vector<Vector2d>>& v, const GraphMetrics& graphMetrics,
    const size_t numX, const size_t numY, const double vScale, const double maxVNorm)
{
    for(size_t i = 0; i < numX; i++)
    {
        for(size_t j = 0; j < numY; j++)
        {
            const auto pixelWidth = (graphMetrics.maxX - graphMetrics.minX) / numX;
            const auto pixelHeight = (graphMetrics.maxY - graphMetrics.minY) / numY;

            const auto arrow = (maxVNorm > 0)
                ? withMaxNorm(vScale * v[i][j], maxVNorm)
                : (vScale * v[i][j]);
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