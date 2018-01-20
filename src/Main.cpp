#include <vector>
#include <math.h>

#include <SDL.h>

struct Vector2f
{
    float x, y;

    Vector2f() {}
    Vector2f(const float x, const float y) : x(x), y(y) {}
};
struct GraphMetrics
{
    float width, height;
    Vector2f pos;
    double minX, maxX;
    double minY, maxY;

    GraphMetrics() {}
};

Vector2f graphPointToPxPoint(const GraphMetrics& graphMetrics, const Vector2f& graphPoint)
{
    const auto xPercentFromLeft = (graphPoint.x - graphMetrics.minX) / (graphMetrics.maxX - graphMetrics.minX);
    const auto yPercentFromBottom = (graphPoint.y - graphMetrics.minY) / (graphMetrics.maxY - graphMetrics.minY);

    return Vector2f(
        graphMetrics.pos.x + (xPercentFromLeft * graphMetrics.width),
        (graphMetrics.pos.y + graphMetrics.height) - (yPercentFromBottom * graphMetrics.height)
    );
}

const int WINDOW_WIDTH = 640;
const int WINDOW_HEIGHT = 480;

void renderLineGraph(SDL_Renderer* renderer, const GraphMetrics& graphMetrics, const double x0, const double dx, const std::vector<double>& ys)
{
    SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);

    for(auto i = 0; i < ys.size() - 1; i++)
    {
        const auto graphPoint0 = Vector2f(x0 + (i * dx), ys[i]);
        const auto pixelPoint0 = graphPointToPxPoint(graphMetrics, graphPoint0);

        const auto graphPoint1 = Vector2f(x0 + ((i + 1) * dx), ys[i + 1]);
        const auto pixelPoint1 = graphPointToPxPoint(graphMetrics, graphPoint1);

        SDL_RenderDrawLine(renderer, (int)pixelPoint0.x, (int)pixelPoint0.y, (int)pixelPoint1.x, (int)pixelPoint1.y);
    }
}

class Step1LinearConvection
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step1LinearConvection()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        // apply initial condition
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            ys[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto c = 0.25f;

        auto newYs = ys;

        for(size_t i = 1; i < ys.size(); i++)
        {
            newYs[i] = ys[i] - (c * (dt / dx) * (ys[i] - ys[i - 1]));
        }

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 41;
    double dx;
    std::vector<double> ys;
};
class Step2NonlinearConvection
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step2NonlinearConvection()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        // apply initial condition
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            ys[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto timeScale = 1.0 / 4;
        auto newYs = ys;

        for(size_t i = 1; i < ys.size(); i++)
        {
            newYs[i] = ys[i] - (ys[i] * ((timeScale * dt) / dx) * (ys[i] - ys[i - 1]));
        }

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 41;
    double dx;
    std::vector<double> ys;
};
class Step3Diffusion
{
public:
    double fixedTimeStep;

    Step3Diffusion()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        const auto sigma = 0.2f;
        fixedTimeStep = sigma * (dx * dx);

        // apply initial condition
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            ys[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto timeScale = 1.0 / 10;
        const auto nu = 0.3f; // viscosity
        auto newYs = ys;

        for(size_t i = 1; i < ys.size() - 1; i++)
        {
            newYs[i] = ys[i] + (nu * ((timeScale * dt) / (dx * dx)) * (ys[i + 1] - (2 * ys[i]) + ys[i - 1]));
        }

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 41;
    double dx;
    std::vector<double> ys;
};
class Step4BurgersEquation
{
public:
    double fixedTimeStep;

    Step4BurgersEquation()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2 * M_PI;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 10;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        ys.resize(numPoints);

        nu = 0.07;
        fixedTimeStep = dx * nu;

        // apply initial condition
        const auto nu = 0.3; // viscosity
        for(size_t i = 0; i < ys.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);

            const auto xMinus2Pi = x - (2 * M_PI);
            const auto e1 = -((x * x) / (4 * nu));
            const auto e2 = -(pow(xMinus2Pi, 2) / (4 * nu));
            const auto phi = exp(e1) + exp(e2);
            const auto dPhiDx = -(1.0 / (2 * nu)) * ((x * exp(e1)) + (xMinus2Pi * exp(e2)));
            ys[i] = -((2 * nu) / phi) * dPhiDx + 4;
        }
    }
    void update(const double dt)
    {
        const auto timeScale = 1.0 / 4;
        const auto scaledDt = timeScale * dt;
        auto newYs = ys;

        for(size_t i = 1; i < ys.size(); i++)
        {
            const auto iPlus1 = (i + 1) % ys.size();

            newYs[i] = ys[i]
                - (((ys[i] * scaledDt) / dx) * (ys[i] - ys[i - 1]))
                + (((nu * scaledDt) / (dx * dx)) * (ys[iPlus1] - (2 * ys[i]) + ys[i - 1]));
        }

        // apply boundary conditions
        newYs[0] = newYs[newYs.size() - 1];

        ys = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, ys);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 101;
    double dx;
    double nu; // viscosity
    std::vector<double> ys;
};

int main(int argc, char* argv[])
{
    bool quit = false;
    
    // init SDL
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window* window = SDL_CreateWindow("12 Steps to Navier-Stokes",
        SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    SDL_Renderer* renderer = SDL_CreateRenderer(window, -1, 0);

    //Step1LinearConvection step;
    //Step2NonlinearConvection step;
    //Step3Diffusion step;
    Step4BurgersEquation step;

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
        while(accumulatedTime >= step.fixedTimeStep)
        {
            step.update(step.fixedTimeStep);
            accumulatedTime -= step.fixedTimeStep;
        }

        // clear window
        SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
        SDL_RenderClear(renderer);

        // draw stored lines
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