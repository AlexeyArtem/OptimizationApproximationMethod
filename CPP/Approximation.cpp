#include "Approximation.h"
#include <stdexcept>
#include <cmath>
#include <vector>
#include "thread"
#include "Matrix.h"


void initArrays(double** points, double** arrArgs, double** arrValues, int sizeArgs, int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex < 0)
		throw std::invalid_argument("");

	for (size_t i = startIndex; i < endIndex; i++)
	{
		for (size_t j = 0; j < sizeArgs; j++)
		{
			arrArgs[i][j] = std::pow(points[i][0], j);
		}
		arrValues[i][0] = points[i][1];
	}
}

void declareArrays(double** arrArgs, double** arrValues, int sizeArgs, int sizeValues, int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex < 0)
		throw std::invalid_argument("");

	for (size_t i = startIndex; i < endIndex; i++)
	{
		arrArgs[i] = new double[sizeArgs];
		arrValues[i] = new double[sizeValues];
	}
}



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
		values[i][0] = point.Y;
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
		values[i][0] = points[i][1];
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

std::list<point>* approximation::parallelMethodOfMinimumRoots(double** points, int countPoints, int degree, double h)
{
	if (degree <= 0)
		throw std::invalid_argument("The degree cannot be a negative number");

	std::list<point>* resultPoints = new std::list<point>();
	int rows = countPoints;

	int argumentsColumns = degree + 1;
	double** arguments = new double* [rows];

	int valuesColumns = 1;
	double** values = new double* [rows];

	for (size_t i = 0; i < rows; i++)
	{
		arguments[i] = new double[argumentsColumns];
		values[i] = new double[valuesColumns];
	}

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	//for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	//{
	//	if (remnant != 0 && i == processor_count - 1)
	//		endIndex += remnant;

	//	threads[i] = std::thread(declareArrays, arguments, values, argumentsColumns, valuesColumns, startIndex, endIndex);

	//	startIndex = endIndex;
	//	endIndex += countElements;
	//}
	//for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	//{
	//	threads[i].join();
	//}

	startIndex = 0;
	endIndex = countElements;
	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(initArrays, points, arguments, values, argumentsColumns, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	/*auto threadsCount = std::thread::hardware_concurrency();
	std::list<std::thread>* threads();*/
	
	/*for (int i = 0; i < rows; i++)
	{
		
		std::thread t1([&]() {
			for (size_t j = 0; j < argumentsColumns; j++)
			{
				arguments[i][j] = std::pow(points[i][0], j);
			}
		});
		std::thread t2([&]() {
			for (size_t j = 0; j < valuesColumns; j++)
			{
				values[i][j] = points[i][1];
			}
		});

		t1.join();
		t2.join();
	}*/
	

	/*for (auto i = threads->begin(); i != threads->end(); i++)
	{
		threads[i]
	}

	for (std::thread thread : threads.begin())
	{
		thread.join();
	}*/

	matrix argumentsMatrix(arguments, rows, argumentsColumns);
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