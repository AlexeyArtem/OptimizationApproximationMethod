#include "matrix.h"
#include <stdexcept>
#include <iostream>
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

matrix matrix::getTransporse() 
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

double matrix::getDeterminant()
{
	if (!isQuadraticity)
		throw new std::exception("Матрица не является квадратичной.");

	matrix triangular = this->getTriangular();
	double det = 1;
	for (int i = 0; i < this->rows; i++)
	{
		det = det * triangular(i, i);
	}

	return det;
}

matrix matrix::getTriangular()
{
	if (!isQuadraticity)
		throw new std::exception("Матрица не является квадратичной.");

	double** values = copyMatrixValues();
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
			if (std::isnan(value))
			{
				values[i][j] = 0;
			}
		}
	}

	return matrix(values, this->rows, this->columns);
}

double** matrix::copyMatrixValues()
{
	double** values = new double* ();
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
	
	return m.getTransporse() * koef;
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

	double algExtra = minor.getDeterminant() * pow(-1, indexRow + indexColumn);

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

matrix matrix::operator+ (matrix matrixB)
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

	//for (size_t i = 0; i < processor_count; i++) // Запуск потоков 
	//{
	//	if (remnant != 0 && i == processor_count - 1)
	//		endIndex += remnant;

	//	threads[i] = std::thread(sum, result, matrixA, matrixB, startIndex, endIndex);

	//	startIndex = endIndex;
	//	endIndex += countElements;
	//}
	//for (size_t i = 0; i < processor_count; i++) // Ожидание завершения потоков 
	//{
	//	threads[i].join();
	//}

	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < this->columns; j++)
		{
			result(i, j) = matrixA(i, j) + matrixB(i, j);
		}
	}
	
	return result;
}

matrix matrix::operator- (matrix matrixB)
{
	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
		throw new std::exception("Матрицы не совпадают по размерности.");

	matrix matrixA = (*this);

	matrix result(this->rows, this->columns);
	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < this->columns; j++)
		{
			result(i, j) = matrixA(i, j) - matrixB(i, j);
		}
	}

	return result;
}

matrix matrix::operator* (matrix matrixB)
{
	if (this->columns != matrixB.rows)
		throw new std::exception("Количество строк и столбцов перемножаемых матриц не совпадают.");

	matrix matrixA = (*this);

	matrix result(matrixA.rows, matrixB.columns);
	for (int i = 0; i < result.rows; i++)
	{
		for (int j = 0; j < result.columns; j++)
		{
			double value = 0;
			for (int k = 0; k < matrixA.columns; k++)
			{
				value += matrixA(i, k) * matrixB(k, j);
			}
			result(i, j) = value;
		}
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