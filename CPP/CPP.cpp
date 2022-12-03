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

	// ����������� ������� �������
	double pivot = arr[(lastIndex + firstIndex) / 2][0];

	while (topOffset <= bottomOffset)
	{
		// � ������� ����� ������� ����������(��������� �� �����) ��������, ������� ������ ������������
		while (arr[topOffset][0] < pivot) {
			topOffset++;
		}
		// � ������ ����� ���������� ��������, ������� ������ ������������
		while (arr[bottomOffset][0] > pivot) {
			bottomOffset--;
		}

		// ������ �������� �������
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
	double** valuesA = new double* [4];
	double** valuesB = new double* [4];
	for (int count = 0; count < 4; ++count)
	{
		valuesA[count] = new double[1, 2, 3, 4];
		valuesB[count] = new double[4, 3, 2, 1];
	}
	
	matrix matrixA = matrix(valuesA, 4, 4);
	double x = 1;
	double& rx = x;
	double y = rx;
	std::cout << rx;
	
	matrix matrixB = matrix(valuesB, 4, 4);
	matrix result = matrixA + matrixB;


	displayMatrix(result, 4, 4);

	int countPoints = 1000000;
	double** points = Generator::generatePoints(1, 30, countPoints);
	quickSortArrRows(points, 0, countPoints - 1);
	//std::list<point>* listPoints = Generator::generatePointsList(1, 9, countPoints);

	int degree = 5;
	double step = 0.001;

	auto start = std::chrono::high_resolution_clock::now();
	std::list<point>* resultApproximation = approximation::methodOfMinimumRoots(points, countPoints, degree, step);
	auto end = std::chrono::high_resolution_clock::now();
	
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);

	//�������� � ������� ������� �� �����.
	std::cout << "time (ms): " << elapsed.count();
}
