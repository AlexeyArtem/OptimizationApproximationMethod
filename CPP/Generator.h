#pragma once
#include <list>
#include "Point.h"
static class Generator
{
public:
	static double** generatePoints(int minValue, int maxValue, int countPoints);
	static std::list<point>* generatePointsList(int minValue, int maxValue, int countPoints);
};

