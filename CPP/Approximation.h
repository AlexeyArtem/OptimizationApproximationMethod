#pragma once
#include <list>
#include "Point.h"
class approximation
{
public:
	static std::list<point>* methodOfMinimumRoots(std::list<point>* points, int degree, double h);
	//static std::list<point>* methodOfMinimumRoots(point* pointsArr, int countPoints, int degree, double h);
};

