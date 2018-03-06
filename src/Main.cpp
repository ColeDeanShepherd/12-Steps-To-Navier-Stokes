#include <vector>
#include <math.h>
#include <sstream>
#include <SDL.h>

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

const auto maxUpdatesPerFrame = 5;

int main(int argc, char* argv[])
{
    bool quit = false;
    
    // init SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("12 Steps to Navier-Stokes",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    //Step1LinearConvection1D step;
    //Step2NonlinearConvection1D step;
    //Step3Diffusion1D step;
    //Step4BurgersEquation1D step;
    //Step5LinearConvection2D step;
    //Step6NonlinearConvection2D step;
    //Step7Diffusion2D step;
    //Step8BurgersEquation2D step;
    //Step9LaplaceEquation2D step;
    //Step10PoissonEquation2D step;
    Step11CavityFlow step;
    //Step12ChannelFlow step;

    auto lastPerfCount = SDL_GetPerformanceCounter();
    double accumulatedTime = 0;

    while(!quit)
    {
        // update timing
        const auto curPerfCount = SDL_GetPerformanceCounter();
        const double dt = (double)(curPerfCount - lastPerfCount) / SDL_GetPerformanceFrequency();
        accumulatedTime += dt;
        lastPerfCount = curPerfCount;

        // handle events
        SDL_Event event;
        while(SDL_PollEvent(&event))
        {
            if(event.type == SDL_QUIT)
            {
                quit = true;
            }
        }

        // update
        auto updatesThisFrame = 0;
        while((updatesThisFrame < maxUpdatesPerFrame) && (accumulatedTime >= step.fixedTimeStep))
        {
            step.update(step.fixedTimeStep);
            accumulatedTime -= step.fixedTimeStep;

            updatesThisFrame++;
        }

        // clear window
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        step.draw(renderer);

        // render window
        SDL_RenderPresent(renderer);
    }

    // cleanup SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}