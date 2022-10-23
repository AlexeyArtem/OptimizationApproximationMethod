#include "matrix.h"
#include <stdexcept>
#include <iostream>
#include <cassert>

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

double& matrix::operator() (int row, int col) 
{
	assert(row >= 0 && row < this->rows);
	assert(col >= 0 && col < this->columns);

	return this->values[row][col];
}

matrix matrix::operator+ (matrix matrixB)
{
	if (matrixB.columns != this->columns || matrixB.rows != this->rows)
		throw new std::exception("Матрицы не совпадают по размерности.");

	matrix matrixA = (*this);

	matrix result(this->rows, this->columns);
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