#pragma once
#include <iomanip>
class matrix
{
private:
	int rows;
	int columns;
	bool isQuadraticity = false;
	double** values;
	void inverse(double** extraMatrix, int startIndex, int endIndex);
	void transpose(matrix* result, int startIndex, int endIndex);
	void determinant(std::atomic<double>& det, matrix* triangular, int startIndex, int endIndex);
	void copy(double** result, int startIndex, int endIndex);
	//void multiplication(matrix* result, matrix* matrixB, int startIndex, int endIndex);

public:
	matrix(double** values, int rows, int columns);
	matrix(int rows, int columns);

	int getRows() { return rows; }
	int getColumns() { return columns; }
	matrix getTranspose();
	matrix getTransposeAsParallel();
	matrix getInverse();
	matrix getInverseAsParallel();
	double getDeterminant();
	double getDeterminantAsParallel();
	matrix getTriangular();
	double getAlgebraicExtra(int indexRow, int indexColumn);
	double** copyMatrixValues();
	double** copyAsParallel();
	void display();
	static matrix multOnMatrix(matrix matrixA, matrix matrixB);
	static matrix multOnNumber(matrix matrixA, double number);

	double& operator() (int i, int j);
	matrix operator+ (matrix matrixB);
	matrix operator- (matrix matrixB);
	matrix operator* (matrix matrixB);
	matrix operator* (double num);
};