#pragma once

#include <SDL.h>

class Step
{
public:
    double fixedTimeStep;

    virtual ~Step() {}
    virtual void update(const double dt) = 0;
    virtual void draw(SDL_Renderer* renderer) = 0;
};