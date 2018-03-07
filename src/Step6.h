#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step6NonlinearConvection2D : public Step
{
public:
    Step6NonlinearConvection2D();
    void applyInitialCondition();
    void applyBoundaryConditions();
    void update(const double dt);
    void draw(SDL_Renderer* renderer);
private:
    const int numX = 81;
    const int numY = 81;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double sigma = 0.2;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> v;
    SDL_Texture* vxHeightMap;
    SDL_Texture* vyHeightMap;
};