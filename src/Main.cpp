#include <vector>

#include <SDL.h>

const int graphWidth = 300;
const int graphHeight = 300;
const int graphX = 10;
const int graphY = 10;
const float graphMinX = 0;
const float graphMaxX = 2;
const float graphMinY = 0;
const float graphMaxY = 2;
const int numPoints = 41;
const float dx = (graphMaxX - graphMinX) / (numPoints - 1);
std::vector<float> ys(numPoints);

void init()
{
    // apply initial condition
    for(auto i = 0; i < ys.size(); i++)
    {
        const auto x = graphMinX + (i * dx);
        ys[i] = ((x >= 0.5f) && (x <= 1)) ? 2 : 1;
    }
}
void update()
{
    const auto dt = 0.0025f;
    const auto c = 1.f;

    auto newYs = ys;

    for(auto i = 1; i < ys.size(); i++)
    {
        newYs[i] = ys[i] - (c * (dt / dx) * (ys[i] - ys[i - 1]));
    }

    ys = newYs;
}
void draw(SDL_Renderer* renderer)
{
    SDL_SetRenderDrawColor(renderer, 128, 128, 128, 255);

    for(auto i = 0; i < numPoints - 1; i++)
    {
        const auto x0 = graphMinX + (i * dx);
        const auto y0 = ys[i];

        const auto px0 = graphX + (((x0 - graphMinX) / (graphMaxX - graphMinX)) * graphWidth);
        const auto py0 = (graphY + graphHeight) - (((y0 - graphMinY) / (graphMaxY - graphMinY)) * graphHeight);

        const auto x1 = graphMinX + ((i + 1) * dx);
        const auto y1 = ys[i + 1];

        const auto px1 = graphX + (((x1 - graphMinX) / (graphMaxX - graphMinX)) * graphWidth);
        const auto py1 = (graphY + graphHeight) - (((y1 - graphMinY) / (graphMaxY - graphMinY)) * graphHeight);

        SDL_RenderDrawLine(renderer, px0, py0, px1, py1);
    }
}

int main(int argc, char * argv[])
{
    bool quit = false;
    
    // init SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("SDL2 line drawing",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 640, 480, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    init();

    while(!quit)
    {
        SDL_Delay(10);

        SDL_Event event;
        SDL_PollEvent(&event);

        if(event.type == SDL_QUIT)
        {
            quit = true;
        }

        // clear window
        SDL_SetRenderDrawColor(renderer, 242, 242, 242, 255);
        SDL_RenderClear(renderer);

        // draw stored lines
        update();
        draw(renderer);

        // render window
        SDL_RenderPresent(renderer);
    }

    // cleanup SDL
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}