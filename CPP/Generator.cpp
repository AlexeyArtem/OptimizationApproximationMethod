#include <stdlib.h>
#include <iostream>
#include "Generator.h"
#include "time.h"

double** Generator::generatePoints(int minValue, int maxValue, int countPoints)
{
	double** points = new double* [countPoints];
	for (int i = 0; i < countPoints; i++)
	{
		points[i] = new double[2];
	}

	srand(time(NULL));
	double randValue;
	int precisionPoints = 2;

	for (int i = 0; i < countPoints; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			randValue = rand() % (int)pow(10, precisionPoints);
			randValue = minValue + (randValue / pow(10, precisionPoints)) * (maxValue - minValue);

			points[i][j] = randValue;
		}
	}

	return points;
}

std::list<point>* Generator::generatePointsList(int minValue, int maxValue, int countPoints)
{
	std::list<point>* points = new std::list<point>();

	srand(time(NULL));
	double value, x, y;
	int precisionPoints = 2;

	for (int i = 0; i < countPoints; i++)
	{
		value = rand() % (int)pow(10, precisionPoints);
		x = minValue + (value / pow(10, precisionPoints)) * (maxValue - minValue);
		value = rand() % (int)pow(10, precisionPoints);
		y = minValue + (value / pow(10, precisionPoints)) * (maxValue - minValue);

		points->push_back(point(x, y));
	}

	return points;
}
