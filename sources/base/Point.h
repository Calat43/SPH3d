#pragma once

class Point {
public:
    double x;
    double y;
    double z;

    Point(double x, double y, double z);

    Point(const Point & p);

    Point operator +(const Point & p);

    Point operator -(const Point & p);

    Point operator *(double term);

    Point operator /(double term);

    Point & operator =(const Point & p);
};
