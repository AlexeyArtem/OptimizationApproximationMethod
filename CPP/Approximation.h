#pragma once
#include <list>
#include "Point.h"
static class approximation
{
public:
	static std::list<point>* methodOfMinimumRoots(std::list<point>* points, int degree, double h);
	static std::list<point>* methodOfMinimumRoots(double** points, int countPoints, int degree, double h);
};

