#pragma once

#include <memory>
#include <SDL.h>
#include "Step.h"

class App
{
public:
    void run();
private:
    const int WINDOW_WIDTH = 640;
    const int WINDOW_HEIGHT = 480;
    const unsigned int MAX_UPDATES_PER_FRAME = 5;

    bool shouldExit;
    Uint64 lastPerfCount;
    double accumulatedTime;
    int stepNumber = 1;
    std::unique_ptr<Step> step;
    SDL_Window* window;
    SDL_Renderer* renderer;

    void initSDL();
    void cleanupSDL();
    std::unique_ptr<Step> createStep(const int stepNumber) const;
    void changeStep(const int newStepNumber);
    void updateTiming();
    void handleSDLEvent(const SDL_Event& event);
    void drawFrame();
    void runFrame();
};