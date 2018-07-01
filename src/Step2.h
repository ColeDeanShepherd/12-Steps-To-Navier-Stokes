#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step2NonlinearConvection1D : public Step
{
public:
    Step2NonlinearConvection1D();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numPoints = 41;

    GraphMetrics graphMetrics;
    double dx;
    std::vector<double> u;
};