using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
using System.Collections;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace CSharp
{
    class Program
    {
        static void Main(string[] args)
        {
            Stopwatch stopWatch = new Stopwatch();

            Matrix matrixA = generateMatrix(1000, 1000, 1, 99);
            Matrix matrixB = generateMatrix(100, 100, 1, 99);
            Matrix matrixC = generateMatrix(10000, 1000, 1, 99);
            Matrix matrixD = generateMatrix(100, 100, 1, 99);
            double number = 33;

            int countRepeats = 3;

            double timeSingleThread, timeManyThread, effectiveness;

            #region TestMNK
                int degree = 5;
                double step = 0.001;
                List<Point> pointsList = Generator.GeneratePoints(1, 30, 1000000);
                //var pointsList = FileService.ReadFromJson<List<Point>>("data.json");

                pointsList = pointsList.AsParallel().OrderBy(p => p.X).ToList();
                //pointsList = pointsList.OrderBy(p => p.X).ToList();
                Point[] pointsArray = pointsList.ToArray();

                // Многопоточный метод
                stopWatch.Start();
                Approximation.ParallelMethodOfMinimumRoots(pointsArray, degree, step);
                stopWatch.Stop();

                timeManyThread = stopWatch.ElapsedMilliseconds;
                Console.WriteLine("Multithreaded MNK time (ms): " + timeManyThread);

                // Однопоточный метод
                stopWatch.Restart();
                Approximation.MethodOfMinimumRoots(pointsArray, degree, step);
                stopWatch.Stop();

                timeSingleThread = stopWatch.ElapsedMilliseconds;
                effectiveness = timeSingleThread / timeManyThread * 100 - 100;
                Console.WriteLine("Single thread MNK time (ms): " + timeSingleThread);
                Console.WriteLine("Effectiveness many thread method = " + effectiveness + "\n");

            #endregion TestMNK

            #region TestDeterminant
            stopWatch.Start();
                for (int i = 0; i < countRepeats; i++)
                {
                    matrixB.GetDeterminantAsParallel();
                }
                stopWatch.Stop();

                timeManyThread = stopWatch.ElapsedMilliseconds;
                Console.WriteLine("Determinant method many thread (ms): " + timeManyThread);

                stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    matrixB.GetDeterminant();
                }
                stopWatch.Stop();

                timeSingleThread = stopWatch.ElapsedMilliseconds;
                effectiveness = timeSingleThread / timeManyThread * 100 - 100;
                Console.WriteLine("Determinant method single thread (ms): " + timeSingleThread);
                Console.WriteLine("Effectiveness many thread method = " + effectiveness + "\n");

            #endregion TestDeterminant

            #region TestTranspose
            stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    matrixC.GetTransporseAsParallel();
                }
                stopWatch.Stop();

                timeManyThread = stopWatch.ElapsedMilliseconds;
                Console.WriteLine("Transpose method many thread (ms): " + timeManyThread);

                stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    matrixC.GetTransporse();
                }
                stopWatch.Stop();

                timeSingleThread = stopWatch.ElapsedMilliseconds;
                effectiveness = timeSingleThread / timeManyThread * 100 - 100;
                Console.WriteLine("Transpose method single thread (ms): " + timeSingleThread);
                Console.WriteLine("Effectiveness many thread method = " + effectiveness + "\n");
            #endregion TestTranspose

            #region TestMultiplicationMatrixOnMatrix
            stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    Matrix.MultiplicationAsParallel(matrixB, matrixD);
                }
                stopWatch.Stop();

                timeManyThread = stopWatch.ElapsedMilliseconds;
                Console.WriteLine("Multiplication method many thread (ms): " + timeManyThread);

                stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    var res = matrixB * matrixD;
                }
                stopWatch.Stop();

                timeSingleThread = stopWatch.ElapsedMilliseconds;
                effectiveness = timeSingleThread / timeManyThread * 100 - 100;
                Console.WriteLine("Multiplication MOM method single thread (ms): " + timeSingleThread);
                Console.WriteLine("Effectiveness many thread method = " + effectiveness + "\n");
            #endregion TestMultiplicationMatrixOnMatrix

            #region TestMultiplicationMatrixOnNumber
            stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    Matrix.MultiplicationAsParallel(matrixB, number);
                }
                stopWatch.Stop();

                timeManyThread = stopWatch.ElapsedMilliseconds;
                Console.WriteLine("Multiplication MON method many thread (ms): " + timeManyThread);

                stopWatch.Restart();
                for (int i = 0; i < countRepeats; i++)
                {
                    var res = matrixB * number;
                }
                stopWatch.Stop();

                timeSingleThread = stopWatch.ElapsedMilliseconds;
                effectiveness = timeSingleThread / timeManyThread * 100 - 100;
                Console.WriteLine("Multiplication MON method single thread (ms): " + timeSingleThread);
                Console.WriteLine("Effectiveness many thread method = " + effectiveness + "\n");
            #endregion TestMultiplicationMatrixOnNumber

            Console.Read();
        }

        static void PrintPoints(IEnumerable<Point> points) 
        {
            foreach (var p in points)
                Console.WriteLine($"X: {p.X}; Y: {p.Y}");
        }

        static Matrix generateMatrix(int rowsCount, int columnsCount, int minValue, int maxValue)
        {
            Random random = new Random();
            double[,] matrix = new double[rowsCount, columnsCount];

            for (int i = 0; i < rowsCount; i++)
            {
                for (int j = 0; j < columnsCount; j++)
                {
                    matrix[i, j] = random.Next(minValue, maxValue) + Math.Round(random.NextDouble(), 2);
                }
            }

            Matrix result = new Matrix(matrix);

            return result;
        }
    }
}
