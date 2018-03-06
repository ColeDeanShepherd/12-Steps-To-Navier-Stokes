#include <vector>
#include <math.h>
#include <sstream>

#include <SDL.h>

namespace std
{
    template <typename T>
    T clamp(T value, T min, T max)
    {
        if(value < min)
        {
            return min;
        }
        else if(value > max)
        {
return max;
        } else
        {
        return value;
        }
    }
}

struct Vector2d
{
    double x, y;

    Vector2d() {}
    Vector2d(const double x, const double y) : x(x), y(y) {}
};

Vector2d negate(const Vector2d& v)
{
    return Vector2d(-v.x, -v.y);
}
Vector2d operator - (const Vector2d& v)
{
    return negate(v);
}

Vector2d add(const Vector2d& a, const Vector2d& b)
{
    return Vector2d(a.x + b.x, a.y + b.y);
}
Vector2d operator + (const Vector2d& a, const Vector2d& b)
{
    return add(a, b);
}

Vector2d sub(const Vector2d& a, const Vector2d& b)
{
    return Vector2d(a.x - b.x, a.y - b.y);
}
Vector2d operator - (const Vector2d& a, const Vector2d& b)
{
    return sub(a, b);
}

double dot(const Vector2d& a, const Vector2d& b)
{
    return (a.x * b.x) + (a.y * b.y);
}
double norm(const Vector2d& v)
{
    return sqrt(dot(v, v));
}
Vector2d div(const Vector2d& v, const double d)
{
    return Vector2d(v.x / d, v.y / d);
}

Vector2d mul(const double s, const Vector2d& v)
{
    return Vector2d(s * v.x, s * v.y);
}
Vector2d operator * (const double s, const Vector2d& v)
{
    return mul(s, v);
}

Vector2d normalized(const Vector2d& v)
{
    const auto vNorm = norm(v);
    return (vNorm > 0) ? div(v, vNorm) : Vector2d(0, 0);
}
Vector2d withMaxNorm(const Vector2d& v, const double maxNorm)
{
    const auto vNorm = norm(v);
    return (vNorm <= maxNorm) ? v : maxNorm * normalized(v);
}

struct GraphMetrics
{
    double width, height;
    Vector2d pos;
    double minX, maxX;
    double minY, maxY;

    GraphMetrics() {}
};

Vector2d graphPointToPxPoint(const GraphMetrics& graphMetrics, const Vector2d& graphPoint)
{
    const auto xPercentFromLeft = (graphPoint.x - graphMetrics.minX) / (graphMetrics.maxX - graphMetrics.minX);
    const auto yPercentFromBottom = (graphPoint.y - graphMetrics.minY) / (graphMetrics.maxY - graphMetrics.minY);

    return Vector2d(
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
        const auto graphPoint0 = Vector2d(x0 + (i * dx), ys[i]);
        const auto pixelPoint0 = graphPointToPxPoint(graphMetrics, graphPoint0);

        const auto graphPoint1 = Vector2d(x0 + ((i + 1) * dx), ys[i + 1]);
        const auto pixelPoint1 = graphPointToPxPoint(graphMetrics, graphPoint1);

        SDL_RenderDrawLine(renderer, (int)pixelPoint0.x, (int)pixelPoint0.y, (int)pixelPoint1.x, (int)pixelPoint1.y);
    }
}

std::vector<double> zeros(const size_t rowCount)
{
    return std::vector<double>(rowCount, 0);
}
std::vector<std::vector<double>> zeros(const size_t rowCount, const size_t columnCount)
{
    return std::vector<std::vector<double>>(rowCount, std::vector<double>(columnCount, 0));
}
std::string toString(const std::vector<std::vector<double>>& matrix)
{
    std::ostringstream sstream;

    for(size_t i = 0; i < matrix.size(); i++)
    {
        if(i > 0)
        {
            sstream << '\n';
        }

        for(size_t j = 0; j < matrix[0].size(); j++)
        {
            if(j > 0)
            {
                sstream << ',';
            }
            
            sstream << matrix[i][j];
        }
    }

    return sstream.str();
}

Vector2d gradient1stOrderCentralDiff(
    const std::vector<std::vector<double>>& f,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto dfdx = (f[xIndex + 1][yIndex] - f[xIndex - 1][yIndex]) / (2 * dx);
    const auto dfdy = (f[xIndex][yIndex + 1] - f[xIndex][yIndex - 1]) / (2 * dy);

    return Vector2d(dfdx, dfdy);
}
double divergence1stOrderBackwardDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto dvxdx = (v[xIndex][yIndex].x - v[xIndex - 1][yIndex].x) / dx;
    const auto dvydy = (v[xIndex][yIndex].y - v[xIndex][yIndex - 1].y) / dy;

    return dvxdx + dvydy;
}
Vector2d laplacian2ndOrderCentralDiff(
    const std::vector<std::vector<Vector2d>>& v,
    const size_t xIndex, const size_t yIndex,
    const double dx, const double dy
)
{
    const auto d2vxdx = (v[xIndex + 1][yIndex].x - (2 * v[xIndex][yIndex].x) + v[xIndex - 1][yIndex].x)
        / (dx * dx);
    const auto d2vxdy = (v[xIndex][yIndex + 1].x - (2 * v[xIndex][yIndex].x) + v[xIndex][yIndex - 1].x)
        / (dy * dy);

    const auto d2vydx = (v[xIndex + 1][yIndex].y - (2 * v[xIndex][yIndex].y) + v[xIndex - 1][yIndex].y)
        / (dx * dx);
    const auto d2vydy = (v[xIndex][yIndex + 1].y - (2 * v[xIndex][yIndex].y) + v[xIndex][yIndex - 1].y)
        / (dy * dy);

    return Vector2d(d2vxdx + d2vxdy, d2vydx + d2vydy);
}

class Step1LinearConvection1D
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step1LinearConvection1D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
        graphMetrics.pos = Vector2d(10, 10);
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
class Step8BurgersEquation2D
{
public:
    double fixedTimeStep;

    Step8BurgersEquation2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        const auto sigma = 0.0009;
        fixedTimeStep = (sigma * dx * dy) / nu;

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
    void update(std::vector<std::vector<double>>& heights, const double dt)
    {
        const auto scaledDt = timeScale * dt;
        const auto& u = uHeights;
        const auto& v = vHeights;
        auto w = heights.data();

        auto newHeights = heights;
        for(size_t i = 1; i < (heights.size() - 1); i++)
        {
            for(size_t j = 1; j < (heights[i].size() - 1); j++)
            {
                newHeights[i][j] = w[i][j]
                    - (((u[i][j] * scaledDt) / dx) * (w[i][j] - w[i - 1][j]))
                    - (((v[i][j] * scaledDt) / dy) * (w[i][j] - w[i][j - 1]))
                    + (((nu * scaledDt) / (dx * dx)) * (w[i + 1][j] - (2 * w[i][j]) + w[i - 1][j]))
                    + (((nu * scaledDt) / (dy * dy)) * (w[i][j + 1] - (2 * w[i][j]) + w[i][j - 1]));
            }
        }

        applyBoundaryConditions(newHeights);
        heights = newHeights;
    }
    void update(const double dt)
    {
        update(uHeights, dt);
        update(vHeights, dt);
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
    const int numX = 41;
    const int numY = 41;
    const double graphMaxHeight = 2;
    const double timeScale = 0.1;
    const double c = 1.0;
    const double nu = 0.01;
    double dx;
    double dy;
    std::vector<std::vector<double>> uHeights;
    std::vector<std::vector<double>> vHeights;
    SDL_Texture* uHeightMap;
    SDL_Texture* vHeightMap;
};
class Step9LaplaceEquation2D
{
public:
    double fixedTimeStep;

    Step9LaplaceEquation2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 1;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        fixedTimeStep = 1.0 / 60;

        initForFunction(p, &heightMap);
    }
    void initForFunction(std::vector<std::vector<double>>& heights, SDL_Texture** heightMap)
    {
        heights = std::vector<std::vector<double>>(numY, std::vector<double>(numX));

        // apply initial condition
        for(size_t rowIndex = 0; rowIndex < heights.size(); rowIndex++)
        {
            for(size_t columnIndex = 0; columnIndex < heights[rowIndex].size(); columnIndex++)
            {
                heights[rowIndex][columnIndex] = 0;
            }
        }

        heightMap = NULL;
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& p)
    {
        const double bCVal = 1;

        // p = 0 at x = 0
        for(size_t i = 0; i < numY; i++)
        {
            p[i][0] = 0;
        }

        // p = y at x = 2
        for(size_t i = 0; i < numY; i++)
        {
            const auto y = graphMetrics.minY + (i * dy);
            p[i][numX - 1] = y;
        }

        // dp/dy = 0 at y = 0
        for(size_t i = 0; i < numX; i++)
        {
            p[0][i] = p[1][i];
        }

        // dp/dy = 0 at y = 1
        for(size_t i = 0; i < numX; i++)
        {
            p[numY - 1][i] = p[numY - 2][i];
        }
    }
    void update(std::vector<std::vector<double>>& p, const double dt)
    {
        const auto dx2 = dx * dx;
        const auto dy2 = dy * dy;

        auto newP = p;
        for(size_t i = 1; i < (p.size() - 1); i++)
        {
            for(size_t j = 1; j < (p[i].size() - 1); j++)
            {
                const auto numerator = (dy2 * (p[i + 1][j] + p[i - 1][j])) + (dx2 * (p[i][j + 1] + p[i][j - 1]));
                const auto denominator = 2 * (dx2 + dy2);
                newP[i][j] = numerator / denominator;
            }
        }

        applyBoundaryConditions(newP);
        p = newP;
    }
    void update(const double dt)
    {
        update(p, dt);
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
        if(heightMap == NULL)
        {
            heightMap = SDL_CreateTexture(renderer, SDL_PIXELFORMAT_RGBA8888, SDL_TEXTUREACCESS_STREAMING, numX, numY);
        }

        updateHeightmap(p, heightMap);

        SDL_Rect uGraphRect;
        uGraphRect.w = (int)graphMetrics.width;
        uGraphRect.h = (int)graphMetrics.height;
        uGraphRect.x = (int)graphMetrics.pos.x;
        uGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, heightMap, NULL, &uGraphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 31;
    const int numY = 31;
    const double graphMaxHeight = 1;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    SDL_Texture* heightMap;
};
class Step10LaplaceEquation2D
{
public:
    double fixedTimeStep;

    Step10LaplaceEquation2D()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 1;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        fixedTimeStep = 1.0 / 60;

        // init p
        p = zeros(numY, numX);

        // init b
        b = zeros(numY, numX);
        b[(int)((double)numY / 4)][(int)((double)numX / 4)] = 100;
        b[(int)(3 * (double)numY / 4)][(int)(3 * (double)numX / 4)] = -100;

        heightMap = NULL;
    }
    void applyBoundaryConditions(std::vector<std::vector<double>>& p)
    {
        const double bCVal = 1;

        // p = 0 at x = 0
        for(size_t i = 0; i < numY; i++)
        {
            p[i][0] = 0;
        }

        // p = 0 at x = 2
        for(size_t i = 0; i < numY; i++)
        {
            p[i][numX - 1] = 0;
        }

        // p = 0 at y = 0
        for(size_t i = 0; i < numY; i++)
        {
            p[0][i] = 0;
        }

        // p = 0 at y = 1
        for(size_t i = 0; i < numY; i++)
        {
            p[numY - 1][i] = 0;
        }
    }
    void update(std::vector<std::vector<double>>& p, const double dt)
    {
        const auto dx2 = dx * dx;
        const auto dy2 = dy * dy;

        auto newP = p;
        for(size_t i = 1; i < (p.size() - 1); i++)
        {
            for(size_t j = 1; j < (p[i].size() - 1); j++)
            {
                const auto numerator = 
                    (dy2 * (p[i + 1][j] + p[i - 1][j])) +
                    (dx2 * (p[i][j + 1] + p[i][j - 1])) - 
                    (b[i][j] * dx2 * dy2);
                const auto denominator = 2 * (dx2 + dy2);
                newP[i][j] = numerator / denominator;
            }
        }

        applyBoundaryConditions(newP);
        p = newP;
    }
    void update(const double dt)
    {
        update(p, dt);
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
                const auto pixelValuePercent =
                    (heights[rowIndex][columnIndex] - graphMinHeight) /
                    (graphMaxHeight - graphMinHeight);
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

        updateHeightmap(p, heightMap);

        SDL_Rect uGraphRect;
        uGraphRect.w = (int)graphMetrics.width;
        uGraphRect.h = (int)graphMetrics.height;
        uGraphRect.x = (int)graphMetrics.pos.x;
        uGraphRect.y = (int)graphMetrics.pos.y;

        SDL_RenderCopy(renderer, heightMap, NULL, &uGraphRect);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 50;
    const int numY = 50;
    const double graphMinHeight = -0.06;
    const double graphMaxHeight = 0.06;
    double dx;
    double dy;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    SDL_Texture* heightMap;
};
class Step11CavityFlow
{
public:
    const double fixedTimeStep = 1.0 / 60.0;

    Step11CavityFlow()
    {
        graphMetrics.width = WINDOW_WIDTH - 20;
        graphMetrics.height = WINDOW_HEIGHT - 20;
        graphMetrics.pos = Vector2d(10, 10);
        graphMetrics.minX = 0;
        graphMetrics.maxX = 2;
        graphMetrics.minY = 0;
        graphMetrics.maxY = 2;

        dx = (graphMetrics.maxX - graphMetrics.minX) / (numX - 1);
        dy = (graphMetrics.maxY - graphMetrics.minY) / (numY - 1);

        b = zeros(numX, numY);
        p = zeros(numX, numY);
        fv = std::vector<std::vector<Vector2d>>(numX, std::vector<Vector2d>(numY, Vector2d(0, 0)));

        heightMap = NULL;
    }

    void buildUpB(const double dt)
    {
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                b[i][j] = getB(i, j, dt);
            }
        }
    }
    double getB(const size_t i, const size_t j, const double dt)
    {
        const auto term1 = (1.0 / dt) * divergence1stOrderBackwardDiff(fv, i, j, dx, dy);
        const auto term2 = pow((fv[i + 1][j].x - fv[i - 1][j].x) / (2 * dx), 2);
        const auto term3 = 2 *
            ((fv[i][j + 1].x - fv[i][j - 1].x) / (2 * dy)) *
            ((fv[i + 1][j].y - fv[i - 1][j].y) / (2 * dx));
        const auto term4 = pow((fv[i][j + 1].y - fv[i][j - 1].y) / (2 * dy), 2);

        return term1 - term2 - term3 - term4;
    }
    void iterateP(const double dt)
    {
        auto newP = p;
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                const auto term1Numerator =
                    ((dy * dy) * (p[i + 1][j] + p[i - 1][j])) +
                    ((dx * dx) * (p[i][j + 1] + p[i][j - 1]));
                const auto termDenominator = 2 * ((dx * dx) + (dy * dy));
                const auto term1 = term1Numerator / termDenominator;
                const auto term2 = ((rho * (dx * dx) * (dy * dy)) / termDenominator) * getB(i, j, dt);

                newP[i][j] = term1 - term2;
            }
        }

        p = newP;
        
        // dp/dy = 0 at y = 0
        for(size_t i = 0; i < numX; i++)
        {
            p[i][0] = p[i][1];
        }

        // p = 0 at y = 2
        for(size_t i = 0; i < numX; i++)
        {
            p[i][numY - 1] = 0;
        }

        // dp/dx = 0 at x = 0, 2
        for(size_t i = 0; i < numX; i++)
        {
            // dp/dx = 0 at x = 0
            p[0][i] = p[1][i];

            // dp/dx = 0 at x = 2
            p[numX - 1][i] = p[numX - 2][i];
        }
    }
    void updateP(const double dt)
    {
        for(size_t iter = 0; iter < numPIterations; iter++)
        {
            iterateP(dt);
        }
    }

    void updateFlowVelocity(const double dt)
    {
        auto newFv = fv;
        for(size_t i = 1; i < numX - 1; i++)
        {
            for(size_t j = 1; j < numY - 1; j++)
            {
                const auto divergenceOfV = divergence1stOrderBackwardDiff(fv, i, j, dx, dy);
                const auto gradientOfP = gradient1stOrderCentralDiff(p, i, j, dx, dy);
                const auto laplacianOfV = laplacian2ndOrderCentralDiff(fv, i, j, dx, dy);

                const auto dfvdt = -(divergenceOfV * fv[i][j])
                    - ((1.0 / rho) * gradientOfP)
                    + (nu * laplacianOfV);

                newFv[i][j] = fv[i][j] + (dt * (dfvdt));
            }
        }

        fv = newFv;

        // v = 0 at boundaries
        for(size_t i = 0; i < numX; i++)
        {
            fv[0][i] = Vector2d(0, 0); // v = 0 at x = 0
            fv[numX - 1][i] = Vector2d(0, 0); // v = 0 at x = 2
            fv[i][0] = Vector2d(0, 0); // v = 0 at y = 0
            fv[i][numY - 1] = Vector2d(0, 0); // v = 0 at y = 2
        }

        // vx = 1 at y = 2
        for(size_t i = 0; i < numX; i++)
        {
            fv[i][numY - 1].x = 1;
        }
    }
    
    void update(const double dt)
    {
        const auto scaledDt = timeScale * dt;

        updateP(scaledDt);
        updateFlowVelocity(scaledDt);
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

        // pixels go from top to bottom, left to right
        for(auto pixelRowIndex = 0; pixelRowIndex < numY; pixelRowIndex++)
        {
            for(auto pixelColumnIndex = 0; pixelColumnIndex < numX; pixelColumnIndex++)
            {
                const auto pixelBytesOffset = (pitch * pixelRowIndex) + (4 * pixelColumnIndex);

                const auto pPercent = std::clamp<double>(
                    (p[pixelColumnIndex][(numY - 1) - pixelRowIndex] - minP) / (maxP - minP),
                    0, 1);

                const auto pixelValue = (unsigned char)(255 * pPercent);
                
                pixelBytes[pixelBytesOffset + 3] = pixelValue; // R
                pixelBytes[pixelBytesOffset + 2] = pixelValue; // G
                pixelBytes[pixelBytesOffset + 1] = pixelValue; // B
                pixelBytes[pixelBytesOffset + 0] = 255; // A
            }
        }

        SDL_UnlockTexture(heightMap);
    }
    void drawVelocityVectors(SDL_Renderer* renderer)
    {
        for(size_t i = 0; i < numX; i++)
        {
            for(size_t j = 0; j < numY; j++)
            {
                const auto pixelWidth = (graphMetrics.maxX - graphMetrics.minX) / numX;
                const auto pixelHeight = (graphMetrics.maxY - graphMetrics.minY) / numY;

                const auto arrow = withMaxNorm(fv[i][j], (3.0 / 4.0) * pixelHeight);
                const auto velocityTailGraphPos = Vector2d(
                    graphMetrics.minX + (i * pixelWidth) + (pixelWidth / 2),
                    graphMetrics.minY + (j * pixelHeight) + (pixelHeight / 2));
                const auto velocityHeadGraphPos = add(velocityTailGraphPos, arrow);
                const auto velocityTailPxPos = graphPointToPxPoint(graphMetrics, velocityTailGraphPos);
                const auto velocityHeadPxPos = graphPointToPxPoint(graphMetrics, velocityHeadGraphPos);

                SDL_SetRenderDrawColor(renderer, 0, 255, 0, 255);
                SDL_RenderDrawLine(
                    renderer, velocityTailPxPos.x, velocityTailPxPos.y,
                    velocityHeadPxPos.x, velocityHeadPxPos.y
                );

                SDL_SetRenderDrawColor(renderer, 255, 0, 0, 255);
                SDL_RenderDrawPoint(renderer, velocityHeadPxPos.x, velocityHeadPxPos.y);
            }
        }
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

        // draw velocity vectors
        drawVelocityVectors(renderer);
    }
private:
    GraphMetrics graphMetrics;
    const int numX = 41;
    const int numY = 41;
    const double timeScale = 0.006;
    const int numPIterations = 50;
    const double c = 1;
    const double rho = 1;
    const double nu = 0.1;
    const double minP = -4;
    const double maxP = 4;
    double dx;
    double dy;
    std::vector<std::vector<Vector2d>> fv;
    std::vector<std::vector<double>> p;
    std::vector<std::vector<double>> b;
    SDL_Texture* heightMap;
};

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
    //Step10LaplaceEquation2D step;
    Step11CavityFlow step;

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
        SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
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