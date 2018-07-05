#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step3Diffusion1D : public Step
{
public:
    Step3Diffusion1D(const int windowWidth, const int windowHeight);
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numPoints = 41;
    const double sigma = 0.2;
    const double nu = 0.3; // viscosity
    const double timeScale = 1.0 / 10;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};