#pragma once

#include <cmath>

struct Vector2d
{
    double x, y;

    Vector2d() {}
    Vector2d(const double x, const double y) : x(x), y(y) {}
};

Vector2d negate(const Vector2d& v);
Vector2d operator - (const Vector2d& v);

Vector2d add(const Vector2d& a, const Vector2d& b);
Vector2d operator + (const Vector2d& a, const Vector2d& b);

Vector2d sub(const Vector2d& a, const Vector2d& b);
Vector2d operator - (const Vector2d& a, const Vector2d& b);

double dot(const Vector2d& a, const Vector2d& b);
double norm(const Vector2d& v);

Vector2d div(const Vector2d& v, const double d);
Vector2d operator / (const Vector2d& v, const double d);

Vector2d mul(const double s, const Vector2d& v);
Vector2d operator * (const double s, const Vector2d& v);

Vector2d normalized(const Vector2d& v);
Vector2d withMaxNorm(const Vector2d& v, const double maxNorm);