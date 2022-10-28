#include "Approximation.h"
#include <stdexcept>
#include <cmath>
#include "Matrix.h"

std::list<point>* approximation::methodOfMinimumRoots(std::list<point>* points, int degree, double h) 
{
	if (degree <= 0)
		throw std::invalid_argument("The degree cannot be a negative number");

	std::list<point>* resultPoints = new std::list<point>();
	int rows = points->size();

	int argColumns = degree + 1;
	double** arguments = new double* [rows];

	int valuesColumns = 1;
	double** values = new double* [rows];

	for (size_t i = 0; i < rows; i++)
	{
		arguments[i] = new double[argColumns];
		values[i] = new double[valuesColumns];
	}

	int i = 0;
	for (point point : *points) 
	{
		for (size_t j = 0; j < argColumns; j++)
		{
			arguments[i][j] = std::pow(point.X, j);
		}
		for (size_t j = 0; j < valuesColumns; j++)
		{
			values[i][j] = point.Y;
		}
		i++;
	}

	matrix argumentsMatrix(arguments, rows, argColumns);
	matrix valuesMatrix(values, rows, valuesColumns);

	matrix first = (argumentsMatrix.getTransporse() * argumentsMatrix).getInverse();
	matrix second = argumentsMatrix.getTransporse() * valuesMatrix;
	matrix A = first * second;

	for (double x = points->front().X; x < points->back().X; x += h)
	{
		double y = 0;
		for (int k = 0; k <= degree; k++)
		{
			y += A(k, 0) * std::pow(x, k);
		}
		resultPoints->push_back(point(x, y));
	}

	// Удалить выделенную память на массивы arguments и values
	/*for (int i = 0; i < sizeY; ++i) {
		delete[] ary[i];
	}
	delete[] ary;*/
	return resultPoints;
}

std::list<point>* approximation::methodOfMinimumRoots(double** points, int countPoints, int degree, double h)
{
	if (degree <= 0)
		throw std::invalid_argument("The degree cannot be a negative number");

	std::list<point>* resultPoints = new std::list<point>();
	int rows = countPoints;

	int argColumns = degree + 1;
	double** arguments = new double* [rows];

	int valuesColumns = 1;
	double** values = new double* [rows];

	for (size_t i = 0; i < rows; i++)
	{
		arguments[i] = new double[argColumns];
		values[i] = new double[valuesColumns];
	}

	for (int i = 0; i < rows; i++)
	{
		for (size_t j = 0; j < argColumns; j++)
		{
			arguments[i][j] = std::pow(points[i][0], j);
		}
		for (size_t j = 0; j < valuesColumns; j++)
		{
			values[i][j] = points[i][1];
		}
	}

	matrix argumentsMatrix(arguments, rows, argColumns);
	matrix valuesMatrix(values, rows, valuesColumns);

	matrix first = (argumentsMatrix.getTransporse() * argumentsMatrix).getInverse();
	matrix second = argumentsMatrix.getTransporse() * valuesMatrix;
	matrix A = first * second;

	for (double x = points[0][0]; x < points[rows - 1][0]; x += h)
	{
		double y = 0;
		for (int k = 0; k <= degree; k++)
		{
			y += A(k, 0) * std::pow(x, k);
		}
		resultPoints->push_back(point(x, y));
	}

	// Удалить выделенную память на массивы arguments и values
	/*for (int i = 0; i < sizeY; ++i) {
		delete[] ary[i];
	}
	delete[] ary;*/
	return resultPoints;
}
