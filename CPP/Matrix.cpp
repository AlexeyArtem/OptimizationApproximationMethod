#include "matrix.h"
#include <stdexcept>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cassert>
#include "thread"

matrix::matrix(double** values, int rows, int columns)
{
	assert(rows > 0 && columns > 0);

	this->values = values;
	this->rows = rows;
	this->columns = columns;
	if (rows == columns) this->isQuadraticity = true;
}

matrix::matrix(int rows, int columns)
{
	assert(rows > 0 && columns > 0);
	this->rows = rows;
	this->columns = columns;
	if (rows == columns) this->isQuadraticity = true;

	double** values = new double* [rows];
	for (int count = 0; count < rows; ++count)
		values[count] = new double[columns];
	this->values = values;
}

matrix matrix::getTranspose() 
{
	matrix result(this->getColumns(), this->getRows());
	matrix matrix = *this;

	for (int i = 0; i < result.getRows(); i++)
	{
		for (int j = 0; j < result.getColumns(); j++)
		{
			result(i, j) = matrix(j, i);
		}
	}

	return result;
}

void matrix::transpose(matrix* result, int startIndex, int endIndex)
{
	for (int i = startIndex; i < endIndex; i++)
	{
		for (int j = 0; j < result->getColumns(); j++)
		{
			(*result)(i, j) = this->values[j][i];
		}
	}
}

matrix matrix::getTransposeAsParallel()
{
	matrix result(this->columns, this->rows);
	
	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = result.rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = result.rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(&matrix::transpose, this, &result, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	return result;
}

void matrix::determinant(std::atomic<double>& det, matrix* triangular, int startIndex, int endIndex)
{
	for (int i = startIndex; i < endIndex; i++)
	{
		det = det * (*triangular)(i, i);
	}
}

double matrix::getDeterminantAsParallel()
{
	if (!isQuadraticity)
		throw new std::exception("Матрица не является квадратичной.");

	matrix triangular = this->getTriangular();

	std::atomic<double> det = 1;

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(&matrix::determinant, this, std::ref(det), &triangular, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	return det;
}

double matrix::getDeterminant()
{
	if (!isQuadraticity)
		throw new std::exception("Матрица не является квадратичной.");

	matrix triangular = this->getTriangular();
	
	double det = 1;
	for (int i = 0; i < this->rows; i++)
	{
		double temp = triangular(i, i);
		det = det * temp;
	}

	return det;
}

matrix matrix::getTriangular()
{
	if (!isQuadraticity)
		throw new std::exception("Матрица не является квадратичной.");

	double** values = copyMatrixValues();

	/*for (size_t i = 0; i < this->rows; i++)
	{
		for (size_t j = 0; j < this->columns; j++)
		{
			std::cout << this->values[i][j] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	for (size_t i = 0; i < this->rows; i++)
	{
		for (size_t j = 0; j < this->columns; j++)
		{
			std::cout << values[i][j] << " ";
		}
		std::cout << "\n";
	}*/

	for (int i = 0; i < this->rows - 1; i++)
	{
		for (int j = i + 1; j < this->columns; j++)
		{
			double coef;
			if (values[i][i] == 0)
			{
				for (int k = 0; k < this->columns; k++)
				{
					values[i][k] += values[i + 1][k];
				}
				coef = values[j][i] / values[i][i];
			}
			else
			{
				coef = values[j][i] / values[i][i];
			}

			//if (std::isnan(coef)) coef = 0;

			for (int k = i; k < this->rows; k++)
			{
				values[j][k] -= values[i][k] * coef;
			}
		}
	}

	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < this->columns; j++)
		{
			double value = values[i][j];
			//std::cout << value << "\n";
			if (std::isnan(value))
			{
				values[i][j] = 0;
			}
		}
	}

	return matrix(values, this->rows, this->columns);
}

void copy(double** result, double** original, int countColumns, int startIndex, int endIndex)
{
	for (size_t i = startIndex; i < endIndex; i++)
	{
		result[i] = new double[countColumns];
		for (size_t j = 0; j < countColumns; j++)
		{
			result[i][j] = original[i][j];
		}
	}
}

double** matrix::copyAsParallel()
{
	double** result = new double* ();

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(copy, result, this->values, this->columns, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	return result;
}

double** matrix::copyMatrixValues()
{
	double** values = new double*[this->rows];
	for (size_t i = 0; i < this->rows; i++)
	{
		values[i] = new double[this->columns];
		for (size_t j = 0; j < this->columns; j++)
		{
			values[i][j] = this->values[i][j];
		}
	}

	return values;
}

void matrix::inverse(double** extraMatrix, int startIndex, int endIndex)
{
	for (int i = startIndex; i < endIndex; i++)
	{
		for (int j = 0; j < this->columns; j++)
		{
			extraMatrix[i][j] = getAlgebraicExtra(i + 1, j + 1);
		}
	}
}

matrix matrix::getInverseAsParallel()
{
	if (!isQuadraticity) throw new std::exception("Матрица не является квадратичной.");

	double det = getDeterminantAsParallel();
	if (det == 0) throw new std::exception("Детерминант равен 0.");

	double** extraMatrix = new double*[this->rows];
	for (size_t i = 0; i < rows; i++)
	{
		extraMatrix[i] = new double[columns];
	}

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(&matrix::inverse, this, extraMatrix, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	double koef = 1 / det;
	matrix m(extraMatrix, rows, columns);

	return m.getTransposeAsParallel() * koef;
}

matrix matrix::getInverse() 
{
	if (!isQuadraticity) throw new std::exception("Матрица не является квадратичной.");

	double det = getDeterminant();
	if (det == 0) throw new std::exception("Детерминант равен 0.");

	double** extraMatrix = new double* [this->rows];
	for (size_t i = 0; i < rows; i++)
	{
		extraMatrix[i] = new double[columns];
	}

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			extraMatrix[i][j] = getAlgebraicExtra(i + 1, j + 1);
		}
	}

	double koef = 1 / det;
	matrix m(extraMatrix, rows, columns);
	
	return m.getTranspose() * koef;
}

double matrix::getAlgebraicExtra(int indexRow, int indexColumn)
{
	if (!isQuadraticity) throw new std::exception("Матрица не является квадратичной.");
	if (indexRow > rows || indexColumn > columns || indexRow < 1 || indexColumn < 1)
		throw new std::exception("Ошибка");

	
	//double** minor = new double* [rows - 1];
	matrix minor = matrix(this->rows - 1, this->columns - 1);

	int row, col;
	row = 0;
	col = 0;

	for (int i = 0; i < this->getRows(); i++)
	{
		for (int j = 0; j < this->getColumns(); j++)
		{
			if (i + 1 != indexRow && j + 1 != indexColumn)
			{
				minor(row, col) = this->values[i][j];
				col++;
			}
		}
		col = 0;
		if (i + 1 != indexRow) row++;
	}

	double algExtra = minor.getDeterminantAsParallel() * pow(-1, indexRow + indexColumn);

	return algExtra;
}

double& matrix::operator() (int row, int col) 
{
	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->columns);

	return this->values[row][col];
}

void sum(matrix result, matrix matrixA, matrix matrixB, int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex < 0)
		throw std::invalid_argument("");

	for (size_t i = startIndex; i < endIndex; i++)
	{
		for (size_t j = 0; j < result.getColumns(); j++)
		{
			result(i, j) = matrixA(i, j) + matrixB(i, j);
		}
	}
}

//matrix matrix::operator+ (matrix matrixB) //Simple threads
//{
//	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
//		throw new std::exception("Матрицы не совпадают по размерности.");
//
//	matrix matrixA = (*this);
//
//	matrix result(this->rows, this->columns);
//
//	for (int i = 0; i < this->rows; i++)
//	{
//		for (int j = 0; j < this->columns; j++)
//		{
//			result(i, j) = matrixA(i, j) + matrixB(i, j);
//		}
//	}
//
//	return result;
//}

matrix matrix::operator+ (matrix matrixB) //Many threads
{
	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
		throw new std::exception("Матрицы не совпадают по размерности.");

	matrix matrixA = (*this);

	matrix result(this->rows, this->columns);

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(sum, result, matrixA, matrixB, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}
	
	return result;
}

void diff(matrix result, matrix matrixA, matrix matrixB, int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex < 0)
		throw std::invalid_argument("");

	for (size_t i = startIndex; i < endIndex; i++)
	{
		for (size_t j = 0; j < result.getColumns(); j++)
		{
			result(i, j) = matrixA(i, j) - matrixB(i, j);
		}
	}
}

//matrix matrix::operator- (matrix matrixB) // Simple thread
//{
//	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
//		throw new std::exception("Матрицы не совпадают по размерности.");
//
//	matrix matrixA = (*this);
//
//	matrix result(this->rows, this->columns);
//	for (int i = 0; i < this->rows; i++)
//	{
//		for (int j = 0; j < this->columns; j++)
//		{
//			result(i, j) = matrixA(i, j) - matrixB(i, j);
//		}
//	}
//
//	return result;
//}

matrix matrix::operator- (matrix matrixB) // Many threads
{
	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
		throw new std::exception("Матрицы не совпадают по размерности.");

	matrix matrixA = (*this);

	matrix result(this->rows, this->columns);

	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(diff, result, matrixA, matrixB, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	return result;
}

//matrix matrix::operator* (matrix matrixB) //Simple thread
//{
//	if (this->columns != matrixB.rows)
//		throw new std::exception("Количество строк и столбцов перемножаемых матриц не совпадают.");
//
//	matrix matrixA = (*this);
//
//	matrix result(matrixA.rows, matrixB.columns);
//	for (int i = 0; i < result.rows; i++)
//	{
//		for (int j = 0; j < result.columns; j++)
//		{
//			double value = 0;
//			for (int k = 0; k < matrixA.columns; k++)
//			{
//				value += matrixA(i, k) * matrixB(k, j);
//			}
//			result(i, j) = value;
//		}
//	}
//
//	return result;
//}

//void matrix::multiplication(matrix* result, matrix* matrixB, int startIndex, int endIndex)
//{
//	if (startIndex < 0 || endIndex < 0)
//		throw std::invalid_argument("");
//
//	matrix matrixA = (*this);
//
//	for (size_t i = startIndex; i < endIndex; i++)
//	{
//		for (size_t j = 0; j < result->getColumns(); j++)
//		{
//			double value = 0;
//			for (int k = 0; k < matrixA.getColumns(); k++)
//			{
//				value += matrixA(i, k) * (*matrixB)(k, j);
//			}
//			(*result)(i, j) = value;
//		}
//	}
//}

void multiplication(matrix result, matrix matrixA, matrix matrixB, int startIndex, int endIndex)
{
	if (startIndex < 0 || endIndex < 0)
		throw std::invalid_argument("");

	for (size_t i = startIndex; i < endIndex; i++)
	{
		for (size_t j = 0; j < result.getColumns(); j++)
		{
			double value = 0;
			for (int k = 0; k < matrixA.getColumns(); k++)
			{
				value += matrixA(i, k) * matrixB(k, j);
			}
			result(i, j) = value;
		}
	}
}

matrix matrix::operator* (matrix matrixB) //Many thread
{
	if (columns != matrixB.rows)
		throw new std::exception("Количество строк и столбцов перемножаемых матриц не совпадают.");

	matrix result(rows, matrixB.columns);
	
	int processor_count = std::thread::hardware_concurrency(); // Количество ядер
	std::thread* threads = new std::thread[processor_count];

	int startIndex = 0;
	int countElements = rows / processor_count; // Количество элементов, которые обрабатываются одним потоком
	int remnant = rows % processor_count; // Остаток от деления
	int endIndex = countElements;

	startIndex = 0;
	endIndex = countElements;

	for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	{
		if (remnant != 0 && i == processor_count - 1)
			endIndex += remnant;

		threads[i] = std::thread(multiplication, result, (*this), matrixB, startIndex, endIndex);

		startIndex = endIndex;
		endIndex += countElements;
	}
	for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	{
		threads[i].join();
	}

	return result;
}

matrix matrix::operator* (double num)
{
	matrix matrixA = (*this);

	matrix result(matrixA.rows, matrixA.columns);
	for (int i = 0; i < result.rows; i++)
	{
		for (int j = 0; j < result.columns; j++)
		{
			result(i, j) = matrixA(i, j) * num;
		}
	}

	return result;
}

void matrix::display()
{
	for (size_t i = 0; i < this->rows; i++)
	{
		for (size_t j = 0; j < this->columns; j++)
		{
			std::cout << std::fixed << std::setprecision(0) << this->values[i][j] << " ";
		}
		std::cout << "\n";
	}
}