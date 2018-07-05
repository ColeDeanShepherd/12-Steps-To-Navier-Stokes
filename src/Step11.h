#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step11CavityFlow : public Step
{
public:
    Step11CavityFlow(const int windowWidth, const int windowHeight);

    void applyPBoundaryConditions();
    void updateP(const double dt);

    void applyFlowVelocityBoundaryConditions();
    void updateFlowVelocity(const double dt);

    void update(const double dt);

    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numX = 41;
    const unsigned int numY = 41;
    const double timeScale = 0.006;
    const unsigned int numPIterations = 50;
    const double c = 1;
    const double rho = 1;
    const double nu = 0.1;
    const double minP = -4;
    const double maxP = 4;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> v;
    std::vector<std::vector<double>> p;
    SDL_Texture* heightMap;
};