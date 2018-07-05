#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step8BurgersEquation2D : public Step
{
public:
    Step8BurgersEquation2D(const int windowWidth, const int windowHeight);
    void applyInitialConditions();
    void applyBoundaryConditions();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double nu = 0.01;
    const double sigma = 0.0009;

    GraphMetrics graphMetrics;
    const unsigned int numX = 41;
    const unsigned int numY = 41;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> v;
    SDL_Texture* vxHeightMap;
    SDL_Texture* vyHeightMap;
};