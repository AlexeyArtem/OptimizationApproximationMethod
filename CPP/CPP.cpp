#include <iostream>
#include <list>
#include <chrono>
#include <time.h> 
#include "matrix.h"
#include "Generator.h"
#include "Approximation.h"
#include "Point.h"

void quickSortArrRows(double** arr, int firstIndex, int lastIndex)
{
	int topOffset = firstIndex;
	int bottomOffset = lastIndex;

	// Центральный элемент массива
	double pivot = arr[(lastIndex + firstIndex) / 2][0];

	while (topOffset <= bottomOffset)
	{
		// В верхней части массива пропускаем(оставляем на месте) элементы, которые меньше центрального
		while (arr[topOffset][0] < pivot) {
			topOffset++;
		}
		// В нижней части пропускаем элементы, которые больше центрального
		while (arr[bottomOffset][0] > pivot) {
			bottomOffset--;
		}

		// Меняем элементы местами
		if (topOffset <= bottomOffset) {
			double* a = arr[topOffset];
			arr[topOffset] = arr[bottomOffset];
			arr[bottomOffset] = a;

			topOffset++;
			bottomOffset--;
		}
	}

	if (bottomOffset > firstIndex) {
		quickSortArrRows(arr, firstIndex, bottomOffset);
	}
	if (topOffset < lastIndex) {
		quickSortArrRows(arr, topOffset, lastIndex);
	}
}

void displayArr(double** arr, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << arr[i][j] << " ";
		}
		std::cout << "\n";
	}
}

void displayMatrix(matrix matrix, int n, int m)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			std::cout << matrix(i, j) << " ";
		}
		std::cout << "\n";
	}
}

void displayList(std::list<point>* points)
{
	for (point point : *points)
	{
		std::cout << "X = " << point.X << "; Y = " << point.Y << "\n";
	}
}

int main()
{
	//double** valuesA = new double* [2];
	//double** valuesB = new double* [2];
	//for (int count = 0; count < 2; count++)
	//{
	//	valuesA[count] = new double[2] {5, 7};
	//	valuesB[count] = new double[2] {3, 5};
	//}

	//matrix matrixA = matrix(valuesA, 2, 2);
	//
	//matrix matrixB = matrix(valuesB, 2, 2);

	//matrix matrixC(3, 3);
	//matrixC(0, 0) = 5;
	//matrixC(0, 1) = 7;
	//matrixC(0, 2) = 3;
	//matrixC(1, 0) = 5;
	//matrixC(1, 1) = 7;
	//matrixC(1, 2) = 2;
	//matrixC(2, 0) = 7;
	//matrixC(2, 1) = 5;
	//matrixC(2, 2) = 3;
	//matrix result = matrixC.getTriangular();

	////std::cout << result;
	//displayMatrix(result, 3, 3);

	int countPoints = 1000000;
	double** points = Generator::generatePoints(1, 30, countPoints);
	quickSortArrRows(points, 0, countPoints - 1);
	//std::list<point>* listPoints = Generator::generatePointsList(1, 9, countPoints);

	int degree = 5;
	double step = 0.001;

	auto start = std::chrono::high_resolution_clock::now();
	std::list<point>* resultApproximation = approximation::parallelMethodOfMinimumRoots(points, countPoints, degree, step);
	auto end = std::chrono::high_resolution_clock::now();
	
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	//Проблема с выводом времени на экран.
	std::cout << "time (ms): " << elapsed.count();
}
