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

matrix generateMatrix(int rowsCount, int columnsCount, int minValue, int maxValue)
{
	matrix resultMatrix(rowsCount, columnsCount);

	srand(time(NULL));
	double randValue;
	int precisionPoints = 2;

	for (size_t i = 0; i < rowsCount; i++)
	{
		for (size_t j = 0; j < columnsCount; j++)
		{
			randValue = rand() % (int)pow(10, precisionPoints);
			randValue = minValue + (randValue / pow(10, precisionPoints)) * (maxValue - minValue);

			resultMatrix(i, j) = randValue;
		}
	}

	return resultMatrix;
}

int main()
{
	int countPoints = 1000000;
	double** points = Generator::generatePoints(1, 30, countPoints);
	quickSortArrRows(points, 0, countPoints - 1);

	int degree = 5;
	double step = 0.001;

	matrix matrixA = generateMatrix(1000, 1000, 1, 99);
	matrix matrixB = generateMatrix(100, 100, 1, 99);
	matrix matrixC = generateMatrix(10, 10, 1, 99);
	matrix matrixD = generateMatrix(100, 100, 1, 99);

	int countRepeats = 50;
	double timeSingleThread, timeManyThread, effectiveness;

	auto start = std::chrono::high_resolution_clock::now();
	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
#pragma region TestTransposeMethods
	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixA.getTranspose();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeSingleThread = elapsed.count() / countRepeats;
	std::cout << "Transpose method single thread (ms): " << timeSingleThread << std::endl;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixA.getTransposeAsParallel();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeManyThread = elapsed.count() / countRepeats;
	effectiveness = timeSingleThread / timeManyThread * 100 - 100;
	std::cout << "Transpose method many thread (ms): " << timeManyThread << std::endl;
	std::cout << "Effectiveness many thread method = " << effectiveness << "%" << std::endl << std::endl;
#pragma endregion

#pragma region TestDeterminantMethods
	countRepeats = 3;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixA.getDeterminant();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeSingleThread = elapsed.count() / countRepeats;
	std::cout << "Determinant method single thread (ms): " << elapsed.count() / countRepeats << std::endl;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixA.getDeterminantAsParallel();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeManyThread = elapsed.count() / countRepeats;
	effectiveness = timeSingleThread / timeManyThread * 100 - 100;
	std::cout << "Determinant method many thread (ms): " << elapsed.count() / countRepeats << std::endl;
	std::cout << "Effectiveness many thread method = " << effectiveness << "%" << std::endl << std::endl;
#pragma endregion

#pragma region TestCopyMethods
	//countRepeats = 50;

	//start = std::chrono::high_resolution_clock::now();
	//for (size_t i = 0; i < countRepeats; i++)
	//{
	//	matrixA.copyMatrixValues();
	//}
	//end = std::chrono::high_resolution_clock::now();

	//elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	//timeSingleThread = elapsed.count() / countRepeats;
	//std::cout << "Copy method single thread (ms): " << elapsed.count() / countRepeats << std::endl;

	//start = std::chrono::high_resolution_clock::now();
	//for (size_t i = 0; i < countRepeats; i++)
	//{
	//	matrixA.copyAsParallel();
	//}
	//end = std::chrono::high_resolution_clock::now();

	//elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	//timeManyThread = elapsed.count() / countRepeats;
	//effectiveness = timeSingleThread / timeManyThread * 100 - 100;
	//std::cout << "Copy method many thread (ms): " << elapsed.count() / countRepeats << std::endl;
	//std::cout << "Effectiveness many thread method = " << effectiveness << "%" << std::endl << std::endl;
#pragma endregion

#pragma region TestInverseMethods
	countRepeats = 10;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixC.getInverse();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeSingleThread = elapsed.count() / countRepeats;
	std::cout << "Inverse method single thread (ms): " << elapsed.count() / countRepeats << std::endl;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		matrixC.getInverseAsParallel();
	}
	end = std::chrono::high_resolution_clock::now();

	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	timeManyThread = elapsed.count() / countRepeats;
	effectiveness = timeSingleThread / timeManyThread * 100 - 100;
	std::cout << "Inverse method many thread (ms): " << elapsed.count() / countRepeats << std::endl;
	std::cout << "Effectiveness many thread method = " << effectiveness << "%" << std::endl << std::endl;
#pragma endregion

#pragma region MNKmethods
	countRepeats = 3;
	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		approximation::methodOfMinimumRoots(points, countPoints, degree, step);
	}
	end = std::chrono::high_resolution_clock::now();

	timeSingleThread = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	std::cout << "Time working single thread method (ms): " << timeSingleThread << std::endl;

	start = std::chrono::high_resolution_clock::now();
	for (size_t i = 0; i < countRepeats; i++)
	{
		approximation::parallelMethodOfMinimumRoots(points, countPoints, degree, step);
	}
	end = std::chrono::high_resolution_clock::now();

	timeManyThread = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
	effectiveness = timeSingleThread / timeManyThread * 100 - 100;
	std::cout << "Time working many threads method (ms): " << timeManyThread << std::endl;
	std::cout << "Effectiveness many thread method = " << effectiveness << "%" << std::endl << std::endl;
#pragma endregion
}
