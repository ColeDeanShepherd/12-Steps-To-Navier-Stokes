#pragma once

#include <string>

struct SDL_Renderer;

class Step
{
public:
    double fixedTimeStep;
    std::string title;

    Step() : title("12 Steps to Navier-Stokes")
    {
    }
    virtual ~Step() {}
    virtual void update(const double dt) = 0;
    virtual void draw(SDL_Renderer* renderer) = 0;
};