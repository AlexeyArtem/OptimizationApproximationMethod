#pragma once
class matrix
{
private:
	int rows;
	int columns;
	bool isQuadraticity = false;
	double** values;

public:
	matrix(double** values, int rows, int columns);
	matrix(int rows, int columns);

	int getRows() { return rows; }
	int getColumns() { return columns; }
	matrix getTransporse();
	matrix getInverse();
	double getDeterminant();
	matrix getTriangular();
	double getAlgebraicExtra(int indexRow, int indexColumn);
	double** copyMatrixValues();

	double& operator() (int i, int j);
	matrix operator+ (matrix matrixB);
	matrix operator- (matrix matrixB);
	matrix operator* (matrix matrixB);
	matrix operator* (double num);
};