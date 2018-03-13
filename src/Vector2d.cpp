#include "Vector2d.h"

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
Vector2d operator / (const Vector2d& v, const double d)
{
    return div(v, d);
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