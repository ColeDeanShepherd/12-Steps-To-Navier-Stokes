#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step4BurgersEquation1D : public Step
{
public:
    Step4BurgersEquation1D();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numPoints = 101;
    const double nu = 0.07; // viscosity
    const double timeScale = 1.0 / 4;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};