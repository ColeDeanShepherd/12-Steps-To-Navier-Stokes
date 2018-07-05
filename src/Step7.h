#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step7Diffusion2D : public Step
{
public:
    Step7Diffusion2D(const int windowWidth, const int windowHeight);
    void applyInitialConditions();
    void applyBoundaryConditions(std::vector<std::vector<double>>& newU);
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numX = 31;
    const unsigned int numY = 31;
    const double graphMaxHeight = 2.0;
    const double timeScale = 0.1;
    const double nu = 0.05;
    const double sigma = 0.25;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> u;
    SDL_Texture* heightMap;
};