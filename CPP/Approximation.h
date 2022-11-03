#pragma once
#include <list>
#include "Point.h"
static class approximation
{
public:
	static std::list<point>* methodOfMinimumRoots(std::list<point>* points, int degree, double h);
	static std::list<point>* methodOfMinimumRoots(double** points, int countPoints, int degree, double h);
	static std::list<point>* parallelMethodOfMinimumRoots(double** points, int countPoints, int degree, double h);
	//void initArrays(std::list<point>* points, double** arrArgs, double** arrValues, int sizeArgs, int startIndex, int endIndex);
	//void declareArrays(double** arrArgs, double** arrValues, int sizeArgs, int sizeValues, int startIndex, int endIndex);
};

