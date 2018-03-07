#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step9LaplaceEquation2D : public Step
{
public:
    Step9LaplaceEquation2D();
    void applyInitialConditions();
    void applyBoundaryConditions(std::vector<std::vector<double>>& p);
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const int numX = 31;
    const int numY = 31;
    const double graphMaxHeight = 1;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    SDL_Texture* heightMap;
};