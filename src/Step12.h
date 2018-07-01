#pragma once

#include <vector>
#include "Step.h"
#include "GraphMetrics.h"

class Step12ChannelFlow : public Step
{
public:
    Step12ChannelFlow();

    void applyPBoundaryConditions();
    void updateP(const double dt);

    void applyFlowVelocityBoundaryConditions();
    void updateFlowVelocity(const double dt);

    void update(const double dt);

    void draw(SDL_Renderer* renderer);
private:
    const unsigned int numX = 41;
    const unsigned int numY = 41;
    const double timeScale = 0.6;
    const unsigned int numPIterations = 50;
    const double c = 1;
    const double rho = 1;
    const double nu = 0.1;
    const double minP = -4;
    const double maxP = 4;

    GraphMetrics graphMetrics;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<Vector2d>> f;
    std::vector<std::vector<Vector2d>> v;
    SDL_Texture* heightMap;
};