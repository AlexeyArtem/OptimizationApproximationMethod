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
            Matrix matrixB = generateMatrix(1000, 1000, 1, 99);
            Matrix matrixC = generateMatrix(1000, 1000, 1, 99);
            Matrix matrixD = generateMatrix(100, 100, 1, 99);
            double number = 33;

            int countRepeats = 3;

            double timeSingleThread, timeManyThread, timeArray, timeList, effectivenessThread, effectivenessStruct;

            #region TestMNK
            Console.WriteLine("The method of least roots:");
            int degree = 5;
            double step = 0.001;
            List<Point> pointsList = Generator.GeneratePoints(1, 30, 1000000);
            //var pointsList = FileService.ReadFromJson<List<Point>>("data.json");

            pointsList = pointsList.AsParallel().OrderBy(p => p.X).ToList();
            Point[] pointsArray = pointsList.ToArray();

            // Однопоточный метод с массивом
            stopWatch.Start();
            Approximation.MethodOfMinimumRoots(pointsArray, degree, step);
            stopWatch.Stop();

            timeArray = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Singlethreaded with array (ms): " + timeArray);

            // Однопоточный метод с листом
            stopWatch.Restart();
            Approximation.MethodOfMinimumRoots(pointsList, degree, step);
            stopWatch.Stop();

            timeList= stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Singlethreaded with list (ms): " + timeList);

            effectivenessStruct = timeList / timeArray;
            Console.WriteLine("Effectiveness struct = " + Math.Round(effectivenessStruct, 3)  + "x\n");
            if (timeArray < timeList) timeSingleThread = timeArray;
            else timeSingleThread = timeList;

            // Многопоточный метод c массивом
            stopWatch.Restart();
            Approximation.ParallelMethodOfMinimumRoots(pointsArray, degree, step);
            stopWatch.Stop();

            timeArray = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Multithreaded with array (ms): " + timeArray);

            // Многопоточный метод c листом
            stopWatch.Restart();
            Approximation.ParallelMethodOfMinimumRoots(pointsList, degree, step);
            stopWatch.Stop();

            timeList = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Multithreaded with list (ms): " + timeList);

            effectivenessStruct = timeList / timeArray;
            Console.WriteLine("Effectiveness struct = " + Math.Round(effectivenessStruct, 3) + "x\n");

            if (timeArray < timeList) timeManyThread = timeArray;
            else timeManyThread = timeList;
            effectivenessThread = timeSingleThread / timeManyThread;
            Console.WriteLine("Effectiveness multithread = " + Math.Round(effectivenessThread, 3) + "x\n");
            #endregion TestMNK

            #region TestDeterminant
            Console.WriteLine("Determinant matrix:");

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                matrixA.GetDeterminant();
            }
            stopWatch.Stop();

            timeSingleThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Singlethreaded (ms): " + timeSingleThread);

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                matrixA.GetDeterminantAsParallel();
            }
            stopWatch.Stop();

            timeManyThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Multithreaded (ms): " + timeManyThread);

            effectivenessThread = timeSingleThread / timeManyThread;
            Console.WriteLine("Effectiveness = " + Math.Round(effectivenessThread, 3) + "x\n");
            #endregion TestDeterminant

            #region TestTranspose
            Console.WriteLine("Transpose matrix:");

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                matrixA.GetTransporse();
            }
            stopWatch.Stop();

            timeSingleThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Singlethreaded (ms): " + timeSingleThread);

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                matrixA.GetTransporseAsParallel();
            }
            stopWatch.Stop();

            timeManyThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Multithreaded (ms): " + timeManyThread);

            effectivenessThread = timeSingleThread / timeManyThread;
            Console.WriteLine("Effectiveness = " + Math.Round(effectivenessThread, 3) + "x\n");
            #endregion TestTranspose

            #region TestMultiplicationMatrixOnMatrix
            Console.WriteLine("Multiplication matrix:");

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                var res = matrixB * matrixC;
            }
            stopWatch.Stop();

            timeSingleThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Singlethreaded (ms): " + timeSingleThread);

            stopWatch.Restart();
            for (int i = 0; i < countRepeats; i++)
            {
                Matrix.MultiplicationAsParallel(matrixB, matrixC);
            }
            stopWatch.Stop();

            timeManyThread = stopWatch.ElapsedMilliseconds;
            Console.WriteLine("Multithreaded (ms): " + timeManyThread);

            effectivenessThread = timeSingleThread / timeManyThread;    
            Console.WriteLine("Effectiveness many thread method = " + Math.Round(effectivenessThread, 3) + "x\n");
            #endregion TestMultiplicationMatrixOnMatrix

            #region TestMultiplicationMatrixOnNumber
            //stopWatch.Restart();
            //    for (int i = 0; i < countRepeats; i++)
            //    {
            //        Matrix.MultiplicationAsParallel(matrixB, number);
            //    }
            //    stopWatch.Stop();

            //    timeManyThread = stopWatch.ElapsedMilliseconds;
            //    Console.WriteLine("Multiplication MON method many thread (ms): " + timeManyThread);

            //    stopWatch.Restart();
            //    for (int i = 0; i < countRepeats; i++)
            //    {
            //        var res = matrixB * number;
            //    }
            //    stopWatch.Stop();

            //    timeSingleThread = stopWatch.ElapsedMilliseconds;
            //    effectivenessThread = timeSingleThread / timeManyThread * 100 - 100;
            //    Console.WriteLine("Multiplication MON method single thread (ms): " + timeSingleThread);
            //    Console.WriteLine("Effectiveness many thread method = " + effectivenessThread + "\n");
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
