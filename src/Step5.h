#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step5LinearConvection2D : public Step
{
public:
    Step5LinearConvection2D();
    void applyInitialConditions();
    void applyBoundaryConditions(std::vector<std::vector<double>>& newHeights);
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double sigma = 0.2;
    const unsigned int numX = 81;
    const unsigned int numY = 81;
    
    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> u;
    SDL_Texture* heightMap;
};