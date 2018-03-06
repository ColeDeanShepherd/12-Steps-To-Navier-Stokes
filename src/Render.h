#pragma once

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