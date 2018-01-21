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

class Step1LinearConvection1D
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step1LinearConvection1D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numPoints - 1);

        heights.resize(numPoints);

        // apply initial condition
        for(size_t i = 0; i < heights.size(); i++)
        {
            const auto x = graphMetrics.minX + (i * dx);
            heights[i] = ((x >= 0.5f) && (x <= 1)) ? 2.0f : 1.0f;
        }
    }
    void update(const double dt)
    {
        const auto c = 0.25f;

        auto newYs = heights;

        for(size_t i = 1; i < heights.size(); i++)
        {
            newYs[i] = heights[i] - (c * (dt / dx) * (heights[i] - heights[i - 1]));
        }

        heights = newYs;
    }
    void draw(SDL_Renderer* renderer)
    {
        renderLineGraph(renderer, graphMetrics, graphMetrics.minX, dx, heights);
    }
private:
    GraphMetrics graphMetrics;
    const int numPoints = 41;
    double dx;
    std::vector<double> heights;
};
class Step2NonlinearConvection1D
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step2NonlinearConvection1D()
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
class Step3Diffusion1D
{
public:
    double fixedTimeStep;

    Step3Diffusion1D()
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
class Step4BurgersEquation1D
{
public:
    double fixedTimeStep;

    Step4BurgersEquation1D()
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
class Step5LinearConvection2D
{
public:
    double fixedTimeStep;

    Step5LinearConvection2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        const auto sigma = 0.2;
        fixedTimeStep = sigma * dx;

        applyInitialConditions();
        heightMap = NULL;
    }
    void applyInitialConditions()
    {
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            const auto y = graphMetrics.minY + (rowIndex * dy);

            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                const auto x = graphMetrics.minX + (columnIndex * dx);
                heights[rowIndex][columnIndex] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
            }
        }
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& newHeights)
    {
        const double bCVal = 1;

        // left
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][0] = bCVal;
        }
        // right
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][numX - 1] = bCVal;
        }
        // top
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[0][i] = bCVal;
        }
        // bottom
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[numY - 1][i] = bCVal;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        auto newHeights = heights;

        // update
        for(size_t i = 1; i < heights.size(); i++)
        {
            for(size_t j = 1; j < heights[i].size(); j++)
            {
                newHeights[i][j] = heights[i][j]
                    - (((c * scaledDt) / dx) * (heights[i][j] - heights[i - 1][j]))
                    - (((c * scaledDt) / dy) * (heights[i][j] - heights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newHeights);

        heights = newHeights;
    }
    void updateHeightmap()
    {
        unsigned char* pixelBytes;
        int pitch;
        const auto failedLocking = SDL_LockTexture(heightMap, NULL, (void**)&pixelBytes, &pitch) != 0;

        if(failedLocking)
        {
            throw new std::exception("Failed locking the texture.");
        }

        for(auto rowIndex = 0; rowIndex < numY; rowIndex++)
        {
            for(auto columnIndex = 0; columnIndex < numX; columnIndex++)
            {
                const auto pixelBytesOffset = (rowIndex * pitch) + (4 * columnIndex);
                const auto pixelValuePercent = heights[rowIndex][columnIndex] / graphMaxHeight;
                const auto pixelValue = (unsigned char)(255 * pixelValuePercent);

                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void draw(SDL_Renderer* renderer)
    {
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap();
        
        SDL_Rect graphRect;
        graphRect.w = (int)graphMetrics.width;
        graphRect.h = (int)graphMetrics.height;
        graphRect.x = (int)graphMetrics.pos.x;
        graphRect.y = (int)graphMetrics.pos.y;
        SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 81;
    const int numY = 81;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    double dx;
    double dy;
    std::vector<std::vector<double>> heights;
    SDL_Texture* heightMap;
};
class Step6NonlinearConvection2D
{
public:
    double fixedTimeStep;

    Step6NonlinearConvection2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        const auto sigma = 0.2;
        fixedTimeStep = sigma * dx;

        initForFunction(uHeights, &uHeightMap);
        initForFunction(vHeights, &vHeightMap);
    }
    void initForFunction(std::vector<std::vector<double>>& heights, SDL_Texture** heightMap)
    {
        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        // apply initial condition
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            const auto y = graphMetrics.minY + (rowIndex * dy);

            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                const auto x = graphMetrics.minX + (columnIndex * dx);
                heights[rowIndex][columnIndex] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
            }
        }

        heightMap = NULL;
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& heights)
    {
        const double bCVal = 1;

        // left
        for(size_t i = 0; i < numY; i++)
        {
            heights[i][0] = bCVal;
        }
        // right
        for(size_t i = 0; i < numY; i++)
        {
            heights[i][numX - 1] = bCVal;
        }
        // top
        for(size_t i = 0; i < numX; i++)
        {
            heights[0][i] = bCVal;
        }
        // bottom
        for(size_t i = 0; i < numX; i++)
        {
            heights[numY - 1][i] = bCVal;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        // update U
        auto newUHeights = uHeights;
        for(size_t i = 1; i < uHeights.size(); i++)
        {
            for(size_t j = 1; j < uHeights[i].size(); j++)
            {
                newUHeights[i][j] = uHeights[i][j]
                    - (((uHeights[i][j] * scaledDt) / dx) * (uHeights[i][j] - uHeights[i - 1][j]))
                    - (((vHeights[i][j] * scaledDt) / dy) * (uHeights[i][j] - uHeights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newUHeights);
        uHeights = newUHeights;

        // update V
        auto newVHeights = vHeights;
        for(size_t i = 1; i < vHeights.size(); i++)
        {
            for(size_t j = 1; j < vHeights[i].size(); j++)
            {
                newVHeights[i][j] = vHeights[i][j]
                    - (((uHeights[i][j] * scaledDt) / dx) * (vHeights[i][j] - vHeights[i - 1][j]))
                    - (((vHeights[i][j] * scaledDt) / dy) * (vHeights[i][j] - vHeights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newVHeights);
        vHeights = newVHeights;
    }
    void updateHeightmap(const std::vector<std::vector<double>>& heights, SDL_Texture* heightMap)
    {
        unsigned char* pixelBytes;
        int pitch;
        const auto failedLocking = SDL_LockTexture(heightMap, NULL, (void**)&pixelBytes, &pitch) != 0;

        if(failedLocking)
        {
            throw new std::exception("Failed locking the texture.");
        }

        for(auto rowIndex = 0; rowIndex < numY; rowIndex++)
        {
            for(auto columnIndex = 0; columnIndex < numX; columnIndex++)
            {
                const auto pixelBytesOffset = (rowIndex * pitch) + (4 * columnIndex);
                const auto pixelValuePercent = heights[rowIndex][columnIndex] / graphMaxHeight;
                const auto pixelValue = (unsigned char)(255 * pixelValuePercent);

                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void draw(SDL_Renderer* renderer)
    {
        if(uHeightMap == NULL)
        {
            uHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        if(vHeightMap == NULL)
        {
            vHeightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap(uHeights, uHeightMap);
        updateHeightmap(vHeights, vHeightMap);

        SDL_Rect uGraphRect;
        uGraphRect.w = ((int)graphMetrics.width / 2) - 5;
        uGraphRect.h = (int)graphMetrics.height;
        uGraphRect.x = (int)graphMetrics.pos.x;
        uGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, uHeightMap, NULL, &uGraphRect);

        SDL_Rect vGraphRect;
        vGraphRect.w = ((int)graphMetrics.width / 2) - 5;
        vGraphRect.h = (int)graphMetrics.height;
        vGraphRect.x = uGraphRect.x + uGraphRect.w + 10;
        vGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, vHeightMap, NULL, &vGraphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 81;
    const int numY = 81;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    double dx;
    double dy;
    std::vector<std::vector<double>> uHeights;
    std::vector<std::vector<double>> vHeights;
    SDL_Texture* uHeightMap;
    SDL_Texture* vHeightMap;
};
class Step7Diffusion2D
{
public:
    double fixedTimeStep;

    Step7Diffusion2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2f(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        const auto sigma = 0.25;
        fixedTimeStep = (sigma * dx * dy) / nu;

        applyInitialConditions();
        heightMap = NULL;
    }
    void applyInitialConditions()
    {
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            const auto y = graphMetrics.minY + (rowIndex * dy);

            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                const auto x = graphMetrics.minX + (columnIndex * dx);
                heights[rowIndex][columnIndex] = ((x >= 0.5) && (x <= 1)) && ((y >= 0.5) && (y <= 1)) ? 2 : 1;
            }
        }
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& newHeights)
    {
        const double bCVal = 1;

        // left
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][0] = bCVal;
        }
        // right
        for(size_t i = 0; i < numY; i++)
        {
            newHeights[i][numX - 1] = bCVal;
        }
        // top
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[0][i] = bCVal;
        }
        // bottom
        for(size_t i = 0; i < numX; i++)
        {
            newHeights[numY - 1][i] = bCVal;
        }
    }
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        auto newHeights = heights;

        // update
        for(size_t i = 1; i < (heights.size() - 1); i++)
        {
            for(size_t j = 1; j < (heights[i].size() - 1); j++)
            {
                newHeights[i][j] = heights[i][j]
                    + (((nu * scaledDt) / (dx * dx)) * (heights[i + 1][j] - (2 * heights[i][j]) + heights[i - 1][j]))
                    + (((nu * scaledDt) / (dy * dy)) * (heights[i][j + 1] - (2 * heights[i][j]) + heights[i][j - 1]));
            }
        }

        applyBoundaryConditions(newHeights);

        heights = newHeights;
    }
    void updateHeightmap()
    {
        unsigned char* pixelBytes;
        int pitch;
        const auto failedLocking = SDL_LockTexture(heightMap, NULL, (void**)&pixelBytes, &pitch) != 0;

        if(failedLocking)
        {
            throw new std::exception("Failed locking the texture.");
        }

        for(auto rowIndex = 0; rowIndex < numY; rowIndex++)
        {
            for(auto columnIndex = 0; columnIndex < numX; columnIndex++)
            {
                const auto pixelBytesOffset = (rowIndex * pitch) + (4 * columnIndex);
                const auto pixelValuePercent = heights[rowIndex][columnIndex] / graphMaxHeight;
                const auto pixelValue = (unsigned char)(255 * pixelValuePercent);

                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void draw(SDL_Renderer* renderer)
    {
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap();

        SDL_Rect graphRect;
        graphRect.w = (int)graphMetrics.width;
        graphRect.h = (int)graphMetrics.height;
        graphRect.x = (int)graphMetrics.pos.x;
        graphRect.y = (int)graphMetrics.pos.y;
        SDL_RenderCopy(renderer, heightMap, NULL, &graphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 31;
    const int numY = 31;
    const double graphMaxHeight = 2.0;
    const double timeScale = 0.1;
    const double nu = 0.05;
    double dx;
    double dy;
    std::vector<std::vector<double>> heights;
    SDL_Texture* heightMap;
};

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
    Step7Diffusion2D step;

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