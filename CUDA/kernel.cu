#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <device_functions.h>

#include <stdio.h>
#include <iostream>
#include <typeinfo>
#include <list>
#include <cmath>
#include <chrono>
#include <time.h> 

using namespace std;



#pragma region Функции ядра
const int THREADS_PER_BLOCK = 1024;

// Функция ядра умножения матриц
__global__ void multMatrixesKernel(double* A, double* B, double* C, int columnsA, int columnsB)
{
    int i0 = columnsA * (blockDim.y * blockIdx.y + threadIdx.y);
    int j0 = blockDim.x * blockIdx.x + threadIdx.x;
    double sum = 0;
    for (int k = 0; k < columnsA; k++)
        sum += A[i0 + k] * B[k * columnsB + j0];
    
    int index = columnsB * (blockDim.y * blockIdx.y + threadIdx.y) + blockDim.x * blockIdx.x + threadIdx.x;
    C[index] = sum;
}

// Функция ядра умножения матрицы на число
__global__ void multMatrixKernel(double* matrix, unsigned int size, double value) 
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < size)
        matrix[index] = matrix[index] * value;
}

// Функция ядра транспонирования матрицы
__global__ void transporseMatrixKernel(double* result, double* matrix, unsigned int rows, unsigned int columns) 
{
    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y
    // (x, y) - координаты, которые уникально идентифицируют поток в сетке

    int i = y * columns + x; // Индекс элемента в исходной матрице
    int j = x * rows + y; // Индекс элемента в транспонированной матрице

    int count = rows * columns;
    if (i < count && j < count) 
    {
        result[j] = matrix[i];
    }
}

__device__ unsigned long long int atomicCAS(unsigned long long int* address, unsigned long long int compare, unsigned long long int val);
__device__ long long int __double_as_longlong(double x);
__device__ double __longlong_as_double(long long int x);
// Атомарная операция умножения числа
__device__ double atomicMult(double* address, double val)
{
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val * __longlong_as_double(assumed)));

    } while (assumed != old);

    return __longlong_as_double(old);
}
__device__ void __syncthreads();

// Функция ядра вычисления детерминанта
__global__ void determinantMatrixKernel(double* triangular, unsigned int size, double* det) 
{
    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size) 
    {
        int index = y * size + y;
        double value = triangular[index];
        atomicMult(det, value);
    }
}

__global__ void reductionDeterminantMatrixKernel_1(double* triangular, unsigned int size, double* resultsBlocks)
{
    __shared__ double blockData[THREADS_PER_BLOCK];

    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size)
    {
        int threadId = threadIdx.x;
        int triangularIndex = y * size + y;

        blockData[threadId] = triangular[triangularIndex];
        __syncthreads();

        for (int i = 1; i < THREADS_PER_BLOCK; i *= 2)
        {
            int index = 2 * i * threadId;
            if ((index + i) < THREADS_PER_BLOCK && blockData[index + i] != 0)
            {
                blockData[index] *= blockData[index + i];
                __syncthreads();
            }
        }
        if (threadId == 0)
            resultsBlocks[blockIdx.x] = blockData[0];
    }
}

__global__ void reductionDeterminantMatrixKernel_2(double* triangular, unsigned int size, double* resultsBlocks)
{
    __shared__ double blockData[THREADS_PER_BLOCK];

    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size)
    {
        int threadId = threadIdx.x;
        int triangularIndex = y * size + y;

        blockData[threadId] = triangular[triangularIndex];
        __syncthreads();

        for (int i = 1; i < THREADS_PER_BLOCK; i *= 2)
        {
            int index = threadId + i;
            if (threadId % (2 * i) == 0 && blockData[index] != 0)
            {
                blockData[threadId] *= blockData[index];
                __syncthreads();
            }
        }
        if (threadId == 0)
            resultsBlocks[blockIdx.x] = blockData[0];
    }
}

__global__ void reductionDeterminantMatrixKernel_3(double* triangular, unsigned int size, double* resultsBlocks)
{
    __shared__ double blockData[THREADS_PER_BLOCK];

    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size)
    {
        int threadId = threadIdx.x;
        int triangularIndex = y * size + y;

        blockData[threadId] = triangular[triangularIndex];
        __syncthreads();

        for (int i = THREADS_PER_BLOCK / 2; i > 0; i >>= 1)
        {
            int index = threadId + i;
            if (threadId < i && blockData[index] != 0)
            {
                blockData[threadId] *= blockData[index];
            }
            __syncthreads();
        }
        if (threadId == 0)
            resultsBlocks[blockIdx.x] = blockData[0];
    }
}

__device__ void warpReduce(volatile double* sdata, int tid)
{
    sdata[tid] += sdata[tid + 32];
    sdata[tid] += sdata[tid + 16];
    sdata[tid] += sdata[tid + 8];
    sdata[tid] += sdata[tid + 4];
    sdata[tid] += sdata[tid + 2];
    sdata[tid] += sdata[tid + 1];
}
__global__ void reductionDeterminantMatrixKernel_5(double* triangular, unsigned int size, double* resultsBlocks)
{
    __shared__ double blockData[THREADS_PER_BLOCK];

    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size)
    {
        int threadId = threadIdx.x;
        int triangularIndex = y * size + y;

        blockData[threadId] = triangular[triangularIndex];
        __syncthreads();

        for (int i = THREADS_PER_BLOCK / 2; i > 0; i >>= 1)
        {
            int index = threadId + i;
            if (threadId < i&& blockData[index] != 0)
            {
                blockData[threadId] *= blockData[index];
            }
            __syncthreads();
        }

        for (unsigned int i = THREADS_PER_BLOCK / 2; i > 32; i >>= 1)
        {
            int index = threadId + i;
            if (threadId < i && blockData[index] != 0)
                blockData[threadId] *= blockData[index];
            __syncthreads();
        }
        if (threadId < 32)
            warpReduce(blockData, threadId);

        if (threadId == 0)
            resultsBlocks[blockIdx.x] = blockData[0];
    }
}

template <unsigned int blockSize>
__device__ void warpReduceTemp(volatile double* sdata, int tid)
{
    if (blockSize >= 64) sdata[tid] += sdata[tid + 32];
    if (blockSize >= 32) sdata[tid] += sdata[tid + 16];
    if (blockSize >= 16) sdata[tid] += sdata[tid + 8];
    if (blockSize >= 8) sdata[tid] += sdata[tid + 4];
    if (blockSize >= 4) sdata[tid] += sdata[tid + 2];
    if (blockSize >= 2) sdata[tid] += sdata[tid + 1];
}
template <unsigned int blockSize>
__global__ void reductionDeterminantMatrixKernel_6(double* triangular, unsigned int size, double* resultsBlocks)
{
    __shared__ double blockData[THREADS_PER_BLOCK];

    int x = blockIdx.x * blockDim.x + threadIdx.x; // Номер потока по X
    int y = blockIdx.y * blockDim.y + threadIdx.y; // Номер потока по Y

    if (y == x && y < size)
    {
        int threadId = threadIdx.x;
        int triangularIndex = y * size + y;

        blockData[threadId] = triangular[triangularIndex];
        __syncthreads();

        for (int i = THREADS_PER_BLOCK / 2; i > 0; i >>= 1)
        {
            int index = threadId + i;
            if (threadId < i&& blockData[index] != 0)
            {
                blockData[threadId] *= blockData[index];
            }
            __syncthreads();
        }

        if (THREADS_PER_BLOCK >= 512) {
            if (threadId < 256) { blockData[threadId] += blockData[threadId + 256]; } __syncthreads();
        }
        if (THREADS_PER_BLOCK >= 256) {
            if (threadId < 128) { blockData[threadId] += blockData[threadId + 128]; } __syncthreads();
        }
        if (THREADS_PER_BLOCK >= 128) {
            if (threadId < 64) { blockData[threadId] += blockData[threadId + 64]; } __syncthreads();
        }
        if (threadId < 32) warpReduceTemp<THREADS_PER_BLOCK>(blockData, threadId);

        if (threadId == 0)
            resultsBlocks[blockIdx.x] = blockData[0];
    }
}


#pragma endregion

void displayDeviceInfo();
void multMatrixesWithCuda(double* C, double* A, double* B, unsigned int rowsA, unsigned int columnsA, unsigned int rowsB, unsigned int columnsB);
void multMatrixWithCuda(double* matrix, unsigned int rows, unsigned int columns, double value);
void transporseMatrixWithCuda(double* result, double* matrix, unsigned int rows, unsigned int columns);
double getDeterminantWithCuda(double* matrix, unsigned int rows, unsigned int columns);
double getDeterminantWithCudaReduction(double* matrix, unsigned int rows, unsigned int columns, int reductionNum, double* elapsedKernel);
double* getTriangularMatrix(double* matrix, unsigned int rows, unsigned int columns);
double* generateMatrix(int rows, int columns, int minValue, int maxValue);

int main()
{
    int countRepeats = 1;
    double time;
    auto start = chrono::high_resolution_clock::now();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double, std::milli> elapsed = end - start;

    unsigned int rowsMatrixA, columnsMatrixA;
    rowsMatrixA = 1000;
    columnsMatrixA = 1000;
    double* matrixA = generateMatrix(rowsMatrixA, columnsMatrixA, 1, 99);

    unsigned int rowsMatrixB, columnsMatrixB;
    rowsMatrixB = 1000;
    columnsMatrixB = 1000;
    double* matrixB = generateMatrix(rowsMatrixB, columnsMatrixB, 1, 99);

    // Тестовый вызов метода для инитиализации среды выполнения CUDA
    cudaSetDevice(0);
    getDeterminantWithCuda(matrixA, rowsMatrixA, columnsMatrixA);
    cout << "===================" << "\n\n";

#pragma region TestTransporseMatrix
    double* resultTransporse = new double[rowsMatrixA * columnsMatrixA];

    start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < countRepeats; i++)
    {
        transporseMatrixWithCuda(resultTransporse, matrixA, rowsMatrixA, columnsMatrixA);
    }
    end = chrono::high_resolution_clock::now();

    elapsed = end - start;
    time = elapsed.count() / countRepeats;
    //cout << "Time transpose method (ms): " << time << "\n\n";
    
    delete[] resultTransporse;
#pragma endregion

#pragma region TestMultMatrixes
    double* resultMultMatrixes = new double[columnsMatrixA * rowsMatrixB];

    start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < countRepeats; i++)
    {
        multMatrixesWithCuda(resultMultMatrixes, matrixA, matrixB, rowsMatrixA, columnsMatrixA, rowsMatrixB, columnsMatrixB);
    }
    end = chrono::high_resolution_clock::now();

    elapsed = end - start;
    time = elapsed.count() / countRepeats;
    //cout << "Time multiplication matrix method (ms): " << time << "\n\n";
    
    delete[] resultMultMatrixes;
#pragma endregion

#pragma region TestDeterminantMatrix
    start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < countRepeats; i++)
    {
        double det = getDeterminantWithCuda(matrixA, rowsMatrixA, columnsMatrixA);
        //cout << "det = " << det << "\n";
    }
    end = chrono::high_resolution_clock::now();

    elapsed = end - start;
    time = elapsed.count() / countRepeats;
    //cout << "Time determinant matrix method (ms): " << time << "\n\n";
#pragma endregion

#pragma region TestReductionDeterminant
    for (int reductionNum = 1; reductionNum <= 6; reductionNum++)
    {
        if (reductionNum == 4) continue;
        double sumElapsedKernel = 0;
        double* elapsedKernel = (double*)malloc(sizeof(double));
        //start = chrono::high_resolution_clock::now();
        for (size_t i = 0; i < countRepeats; i++)
        {
            getDeterminantWithCudaReduction(matrixA, rowsMatrixA, columnsMatrixA, reductionNum, elapsedKernel);
            sumElapsedKernel += *elapsedKernel;
        }
        //end = chrono::high_resolution_clock::now();

        //elapsed = end - start;
        //time = elapsed.count() / countRepeats;
        time = sumElapsedKernel / countRepeats;
        cout << "Time reduction " << reductionNum << " determinant matrix kernel method(ms): " << time << "\n\n";
    }

    /*elapsed = end - start;
    time = elapsed.count() / countRepeats;
    cout << "Time reduction determinant matrix method (ms): " << time << "\n\n";*/
#pragma endregion

#pragma region TestMultMatrix
    start = chrono::high_resolution_clock::now();
    for (size_t i = 0; i < countRepeats; i++)
    {
        multMatrixWithCuda(matrixA, rowsMatrixA, columnsMatrixA, 2);
    }
    end = chrono::high_resolution_clock::now();

    elapsed = end - start;
    time = elapsed.count() / countRepeats;
    //cout << "Time multiplication matrix by number method (ms): " << time << "\n\n";
#pragma endregion

    return 0;
}

// Вывод информации об устройстве (GPU)
void displayDeviceInfo()
{
    const int kb = 1024;
    const int mb = kb * kb;
    std::wcout << "NBody.GPU" << std::endl << "=========" << std::endl << std::endl;

    wcout << "CUDA version:   v" << CUDART_VERSION << endl;
    //wcout << "Thrust version: v" << THRUST_MAJOR_VERSION << "." << THRUST_MINOR_VERSION << endl << endl;

    int devCount;
    cudaGetDeviceCount(&devCount);
    wcout << "CUDA Devices: " << endl << endl;

    for (int i = 0; i < devCount; ++i)
    {
        cudaDeviceProp props;
        cudaGetDeviceProperties(&props, i);
        wcout << i << ": " << props.name << ": " << props.major << "." << props.minor << endl;
        wcout << "  Global memory:   " << props.totalGlobalMem / mb << "mb" << endl;
        wcout << "  Shared memory:   " << props.sharedMemPerBlock / kb << "kb" << endl;
        wcout << "  Constant memory: " << props.totalConstMem / kb << "kb" << endl;
        wcout << "  Block registers: " << props.regsPerBlock << endl << endl;

        wcout << "  Warp size:         " << props.warpSize << endl;
        wcout << "  Threads per block: " << props.maxThreadsPerBlock << endl;
        wcout << "  Max block dimensions: [ " << props.maxThreadsDim[0] << ", " << props.maxThreadsDim[1] << ", " << props.maxThreadsDim[2] << " ]" << endl;
        wcout << "  Max grid dimensions:  [ " << props.maxGridSize[0] << ", " << props.maxGridSize[1] << ", " << props.maxGridSize[2] << " ]" << endl;
        wcout << endl;
    }
}

// Транспонирование матрицы
void transporseMatrixWithCuda(double* result, double* matrix, unsigned int rows, unsigned int columns)
{
    // Выделение памяти на устройстве
    double* devMatrix;
    double* devResult;
    cudaMalloc((void**)&devMatrix, rows * columns * sizeof(double));
    cudaMalloc((void**)&devResult, rows * columns * sizeof(double));

    // Копирование значений с памяти хоста на устройство
    cudaMemcpy(devMatrix, matrix, rows * columns * sizeof(double), cudaMemcpyHostToDevice);

    // Вычисление размеров блока и сетки
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 0);
    int threadsPerBlockDim = sqrt(props.maxThreadsPerBlock); // Max threads for block = 1024 = 32*32
    dim3 blockDim(threadsPerBlockDim, threadsPerBlockDim); // 2D-block 32*32

    int blocksPerGridDimX = ceilf(columns / (float)threadsPerBlockDim);
    int blocksPerGridDimY = ceilf(rows / (float)threadsPerBlockDim);
    dim3 gridDim(blocksPerGridDimX, blocksPerGridDimY);

    auto start = chrono::high_resolution_clock::now();

    // Запуск функции ядра
    transporseMatrixKernel<<< gridDim, blockDim >>>(devResult, devMatrix, rows, columns);
    cudaThreadSynchronize();

    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    cout << "Time transporse matrix kernel method (ms): " << elapsed.count() << std::endl;

    cudaMemcpy(result, devResult, rows * columns * sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(devResult);
    cudaFree(devMatrix);
}

// Умножение матриц
void multMatrixesWithCuda(double* C, double* A, double* B, unsigned int rowsA, unsigned int columnsA, unsigned int rowsB, unsigned int columnsB)
{
    if (columnsA != rowsB)
        throw new std::exception("Количество строк и столбцов перемножаемых матриц не совпадают.");

    double *devC, *devA, *devB;

    // Выделение памяти на устройстве
    cudaMalloc((void**)&devC, columnsA * rowsB * sizeof(double));
    cudaMalloc((void**)&devA, columnsA * rowsA * sizeof(double));
    cudaMalloc((void**)&devB, columnsB * rowsB * sizeof(double));

    // Копирование данных с хоста на устройство
    cudaMemcpy(devA, A, columnsA * rowsA * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devB, B, columnsB * rowsB * sizeof(double), cudaMemcpyHostToDevice);
    
    // Вычисление размеров блока и сетки
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 0);
    int threadsPerBlockDim = sqrt(props.maxThreadsPerBlock);
    dim3 blockDim(threadsPerBlockDim, threadsPerBlockDim);

    int blocksPerGridDimX = ceilf(columnsA / threadsPerBlockDim);
    int blocksPerGridDimY = ceilf(rowsB / threadsPerBlockDim);
    dim3 gridDim(blocksPerGridDimX, blocksPerGridDimY);

    auto start = chrono::high_resolution_clock::now();
    
    // Запуск функции ядра
    multMatrixesKernel << < gridDim, blockDim >> > (devA, devB, devC, columnsA, columnsB);
    cudaThreadSynchronize();

    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    cout << "Time multiplication matrixes kernel method (ms): " << elapsed.count() << std::endl;

    // Копирование результата с устройства на хост
    cudaMemcpy(C, devC, columnsA * rowsB * sizeof(double), cudaMemcpyDeviceToHost);

    // Освобождение памяти устройства
    cudaFree(devA);
    cudaFree(devB);
    cudaFree(devC);
}

// Умножение матрицы на число
void multMatrixWithCuda(double* matrix, unsigned int rows, unsigned int columns, double value) 
{
    int size = rows * columns;

    double* devMatrix;
    cudaMalloc((void**)&devMatrix, size * sizeof(double));
    cudaMemcpy(devMatrix, matrix, size * sizeof(double), cudaMemcpyHostToDevice);

    // Вычисление размеров блока и сетки
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 0);
    int threadsPerBlock = props.maxThreadsPerBlock;
    int blocksPerGrid = (size + threadsPerBlock - 1) / threadsPerBlock;

    auto start = chrono::high_resolution_clock::now();

    // Запуск функции ядра
    multMatrixKernel<<<blocksPerGrid, threadsPerBlock>>>(devMatrix, size, value);
    cudaThreadSynchronize();

    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    cout << "Time multiplication matrix by number kernel method (ms): " << elapsed.count() << std::endl;

    // Копирование результата с устройства на хост
    cudaMemcpy(matrix, devMatrix, size * sizeof(double), cudaMemcpyDeviceToHost);

    // Освобождение памяти устройства
    cudaFree(devMatrix);
}

// Получение треугольной матрицы
double* getTriangularMatrix(double* matrix, unsigned int rows, unsigned int columns) 
{
    if (rows != columns)
        throw std::invalid_argument("Матрица не является квадратичной.");

    double* result = new double[rows * columns];
    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            result[i * columns + j] = matrix[i * columns + j];
        }
    }

    for (int i = 0; i < rows - 1; i++)
    {
        for (int j = i + 1; j < columns; j++)
        {
            double coef;
            if (result[i * columns + i] == 0)
            {
                for (int k = 0; k < columns; k++)
                {
                    result[i * columns + k] += result[(i + 1) * columns + k];
                }
                coef = result[j * columns + i] / result[i * columns + i];
            }
            else
            {
                coef = result[j * columns + i] / result[i * columns + i];
            }

            if (std::isnan(coef)) 
                coef = 0;

            for (int k = i; k < rows; k++)
            {
                result[j * rows + k] -= result[i * rows + k] * coef;
            }
        }
    }

    return result;
}

// Получение детерминанта
double getDeterminantWithCuda(double* matrix, unsigned int rows, unsigned int columns)
{
    if (rows != columns) 
        throw std::invalid_argument("Матрица не является квадратичной.");

    double* triangular = getTriangularMatrix(matrix, rows, columns);
    double* det = (double*)malloc(sizeof(double));
    *det = 1;

    // Выделение памяти на устройстве
    double* devTriangular;
    double* devDet;
    cudaMalloc((void**)&devTriangular, columns * rows * sizeof(double));
    cudaMalloc((void**)&devDet, sizeof(double));
    
    cudaMemcpy(devTriangular, triangular, rows * columns * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(devDet, det, sizeof(double), cudaMemcpyHostToDevice);

    // Вычисление размеров блока и сетки
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 0);
    int threadsPerBlockDim = sqrt(props.maxThreadsPerBlock);
    dim3 blockDim(threadsPerBlockDim, threadsPerBlockDim);

    int blocksPerGridDimX = ceilf(rows / (float)threadsPerBlockDim);
    int blocksPerGridDimY = ceilf(rows / (float)threadsPerBlockDim);
    dim3 gridDim(blocksPerGridDimX, blocksPerGridDimY);

    auto start = chrono::high_resolution_clock::now();

    // Запуск ядра
    determinantMatrixKernel<<<gridDim, blockDim>>> (devTriangular, rows, devDet);
    cudaThreadSynchronize();

    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    cout << "Time determinant matrix kernel method (ms): " << elapsed.count() << std::endl;

    cudaMemcpy(det, devDet, sizeof(double), cudaMemcpyDeviceToHost);

    // Освобождение памяти на хосте и на устройстве
    double result = *det;

    cudaFree(devTriangular);
    cudaFree(devDet);
    
    delete[] triangular;
    free(det);
    
    return result;
}

// Получение детерминанта с помощью редукции
double getDeterminantWithCudaReduction(double* matrix, unsigned int rows, unsigned int columns, int reductionNum, double* elapsedKernel)
{
    if (rows != columns)
        throw std::invalid_argument("Матрица не является квадратичной.");

    double* triangular = getTriangularMatrix(matrix, rows, columns);

    int size = rows;
    cudaDeviceProp props;
    cudaGetDeviceProperties(&props, 0);
    int threadsPerBlock = props.maxThreadsPerBlock;
    //int blocks = ceilf((size + threadsPerBlock - 1) / (float)threadsPerBlock);

    int threadsPerBlockDim = sqrt(threadsPerBlock);
    dim3 blockDim(threadsPerBlockDim, threadsPerBlockDim);

    int blocksPerGridDim = ceilf(rows / (float)threadsPerBlockDim);
    int blocks = blocksPerGridDim * blocksPerGridDim;
    dim3 gridDim(blocksPerGridDim, blocksPerGridDim);

    // Выделение памяти на устройстве
    double* devTriangular;
    double* devResultsBlocks; // Результаты вычисления блоков на устройстве
    cudaMalloc((void**)&devTriangular, rows * columns * sizeof(double));
    cudaMalloc((void**)&devResultsBlocks, blocks * sizeof(double));
    cudaMemcpy(devTriangular, triangular, rows * columns * sizeof(double), cudaMemcpyHostToDevice);

    auto start = chrono::high_resolution_clock::now();
    // Запуск ядра
    switch (reductionNum)
    {
        case 1:
            reductionDeterminantMatrixKernel_1 << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);

        case 2:
            reductionDeterminantMatrixKernel_2 << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);
        case 3:
            reductionDeterminantMatrixKernel_3 << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);
        case 5:
            reductionDeterminantMatrixKernel_5 << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);
        case 6:
            reductionDeterminantMatrixKernel_6<THREADS_PER_BLOCK> << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);

        default:
            reductionDeterminantMatrixKernel_1 << <gridDim, blockDim >> > (devTriangular, size, devResultsBlocks);
    }
    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    *elapsedKernel = elapsed.count();

    cudaDeviceSynchronize();

    double* resultsBlocks = (double*)malloc(blocks * sizeof(double)); // Результаты вычисления блоков на хосте
    cudaMemcpy(resultsBlocks, devResultsBlocks, blocks * sizeof(double), cudaMemcpyDeviceToHost);

    // Вычисление результатов блоков
    double det = 1;
    for (size_t i = 0; i < blocks; i++)
    {
        if (resultsBlocks[i] != 0)
            det *= resultsBlocks[i];
    }

    // Освобождение памяти
    cudaFree(devTriangular);
    cudaFree(devResultsBlocks);
    delete[] triangular;
    free(resultsBlocks);

    return det;
}

// Генерации матрицы случайными числами
double* generateMatrix(int rows, int columns, int minValue, int maxValue) 
{
    double* matrix = new double[rows * columns];
    srand(time(NULL));
    double randValue;
    int precisionPoints = 2;

    for (size_t i = 0; i < rows; i++)
    {
        for (size_t j = 0; j < columns; j++)
        {
            randValue = rand() % (int)pow(10, precisionPoints);
            randValue = minValue + (randValue / pow(10, precisionPoints)) * (maxValue - minValue);
            matrix[i * columns + j] = randValue;
        }
    }

    return matrix;
}