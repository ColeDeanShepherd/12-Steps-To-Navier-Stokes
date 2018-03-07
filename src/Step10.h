#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step10PoissonEquation2D : public Step
{
public:
    Step10PoissonEquation2D();
    void applyBoundaryConditions();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const int numX = 50;
    const int numY = 50;
    const double graphMinHeight = -0.06;
    const double graphMaxHeight = 0.06;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    SDL_Texture* heightMap;
};