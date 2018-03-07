#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step1LinearConvection1D : public Step
{
public:
    Step1LinearConvection1D();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const int numPoints = 41;
    const double c = 0.25;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};