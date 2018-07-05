#include "App.h"
#include "Core.h"
#include "Step1.h"
#include "Step2.h"
#include "Step3.h"
#include "Step4.h"
#include "Step5.h"
#include "Step6.h"
#include "Step7.h"
#include "Step8.h"
#include "Step9.h"
#include "Step10.h"
#include "Step11.h"
#include "Step12.h"

void App::run()
{
    initSDL();
    changeStep(1);

    lastPerfCount = SDL_GetPerformanceCounter();
    accumulatedTime = 0;

#ifndef __EMSCRIPTEN__
    shouldExit = false;
    while(!shouldExit)
    {
        runFrame();
    }
#else
    emscripten_set_main_loop(runFrame, 0, 1);
#endif

    cleanupSDL();
}
void App::initSDL()
{
    SDL_Init(SDL_INIT_VIDEO);
    window = SDL_CreateWindow(
        "12 Steps to Navier-Stokes", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
        WINDOW_WIDTH, WINDOW_HEIGHT, 0
    );
    renderer = SDL_CreateRenderer(window, -1, 0);
}
void App::cleanupSDL()
{
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
}
std::unique_ptr<Step> App::createStep(const int stepNumber) const
{
    switch(stepNumber)
    {
    case 1:
        return std::make_unique<Step1LinearConvection1D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 2:
        return std::make_unique<Step2NonlinearConvection1D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 3:
        return std::make_unique<Step3Diffusion1D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 4:
        return std::make_unique<Step4BurgersEquation1D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 5:
        return std::make_unique<Step5LinearConvection2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 6:
        return std::make_unique<Step6NonlinearConvection2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 7:
        return std::make_unique<Step7Diffusion2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 8:
        return std::make_unique<Step8BurgersEquation2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 9:
        return std::make_unique<Step9LaplaceEquation2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 10:
        return std::make_unique<Step10PoissonEquation2D>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 11:
        return std::make_unique<Step11CavityFlow>(WINDOW_WIDTH, WINDOW_HEIGHT);
    case 12:
        return std::make_unique<Step12ChannelFlow>(WINDOW_WIDTH, WINDOW_HEIGHT);
    default:
        throw std::runtime_error("Invalid step number.");
    }
}
void App::changeStep(const int newStepNumber)
{
    stepNumber = newStepNumber;
    step = createStep(stepNumber);
    SDL_SetWindowTitle(window, step->title.c_str());
}
void App::updateTiming()
{
    const auto curPerfCount = SDL_GetPerformanceCounter();
    const double dt = (double)(curPerfCount - lastPerfCount) / SDL_GetPerformanceFrequency();
    accumulatedTime += dt;
    lastPerfCount = curPerfCount;
}
void App::handleSDLEvent(const SDL_Event& event)
{
    if(event.type == SDL_QUIT)
    {
        shouldExit = true;
    }
    else if(event.type == SDL_KEYDOWN)
    {
        if(event.key.keysym.sym == SDLK_LEFT)
        {
            changeStep(wrap(stepNumber - 1, 1, 12));
        }
        else if(event.key.keysym.sym == SDLK_RIGHT)
        {
            changeStep(wrap(stepNumber + 1, 1, 12));
        }
    }
}
void App::drawFrame()
{
    SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
    SDL_RenderClear(renderer);

    step->draw(renderer);

    SDL_RenderPresent(renderer);
}
void App::runFrame()
{
    updateTiming();

    SDL_Event event;
    while(SDL_PollEvent(&event))
    {
        handleSDLEvent(event);
    }

    auto updatesThisFrame = 0;
    while((updatesThisFrame < MAX_UPDATES_PER_FRAME) && (accumulatedTime >= step->fixedTimeStep))
    {
        step->update(step->fixedTimeStep);
        accumulatedTime -= step->fixedTimeStep;

        updatesThisFrame++;
    }

    drawFrame();
}