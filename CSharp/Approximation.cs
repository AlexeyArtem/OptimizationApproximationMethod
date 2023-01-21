using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;

namespace CSharp
{
    static class Approximation
    {
        public static List<Point> MethodOfMinimumRoots(List<Point> points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Count, degree + 1]; // X
            double[,] values = new double[points.Count, 1]; // Y

            int i = 0;
            foreach (var point in points)
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(point.X, j);
                }
                values[i, 0] = point.Y;
                i++;
            }

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = (argumentsMatrix.GetTransporse() * argumentsMatrix).GetInverse();
            Matrix second = argumentsMatrix.GetTransporse() * valuesMatrix;
            Matrix A = first * second;

            for (double x = points.First().X; x <= points.Last().X; x += h)
            {
                double y = 0;
                for (int k = 0; k <= degree; k++)
                {
                    y += A[k, 0] * Math.Pow(x, k);
                }
                resultPoints.Add(new Point(x, y));
            }

            return resultPoints;
        }

        public static List<Point> MethodOfMinimumRoots(Point[] points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Length, degree + 1]; // X
            double[,] values = new double[points.Length, 1]; // Y

            for (int i = 0; i < points.Length; i++)
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(points[i].X, j);
                }
                values[i, 0] = points[i].Y;
            }

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = (argumentsMatrix.GetTransporse() * argumentsMatrix).GetInverse();
            Matrix second = argumentsMatrix.GetTransporse() * valuesMatrix;
            Matrix A = first * second;

            for (double x = points[0].X; x <= points[points.Length - 1].X; x += h)
            {
                double y = 0;
                for (int k = 0; k <= degree; k++)
                {
                    y += A[k, 0] * Math.Pow(x, k);
                }
                resultPoints.Add(new Point(x, y));
            }

            return resultPoints;
        }

        public static List<Point> ParallelMethodOfMinimumRoots(List<Point> points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Count, degree + 1]; // X
            double[,] values = new double[points.Count, 1]; // Y

            #region Старая реализация потоков
            //int countRowsArgs = arguments.GetLength(0);
            //int processorCount = Environment.ProcessorCount;

            //Initialization[] initsArgs = new Initialization[processorCount % 2 == 0 ? (processorCount / 2) : (processorCount / 2 + 1)];
            //Thread[] threadsArgs = new Thread[initsArgs.Length];

            //int startIndex = 0;
            //int endIndex = countRowsArgs / initsArgs.Length;
            //for (int i = 0; i < initsArgs.Length; i++)
            //{
            //    initsArgs[i] = new Initialization(arguments, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRowsArgs / initsArgs.Length;
            //    if (endIndex == countRowsArgs - 1)
            //        endIndex += 1;

            //    threadsArgs[i] = new Thread(initsArgs[i].InitializeArguments);
            //    threadsArgs[i].Start();
            //}

            //int countRowsValues = values.GetLength(0);
            //Initialization[] initsValues = new Initialization[processorCount / 2];
            //Thread[] threadsValues = new Thread[initsValues.Length];

            //startIndex = 0;
            //endIndex = countRowsValues / initsValues.Length;

            //for (int i = 0; i < initsValues.Length; i++)
            //{
            //    initsValues[i] = new Initialization(values, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRowsValues / initsValues.Length;
            //    if (endIndex == countRowsArgs - 1)
            //        endIndex += 1;

            //    threadsValues[i] = new Thread(initsValues[i].InitializeValues);
            //    threadsValues[i].Start();
            //}

            //foreach (var thread in threadsArgs.Concat(threadsValues))
            //{
            //    thread.Join();
            //}
            #endregion

            #region Новая реализация потоков
            //Initialization[] inits = new Initialization[Environment.ProcessorCount];
            //Thread[] threads = new Thread[Environment.ProcessorCount];

            //int countRows = points.Count;
            //int startIndex = 0;
            //int endIndex = countRows / inits.Length;
            //for (int i = 0; i < inits.Length; i++)
            //{
            //    inits[i] = new Initialization(arguments, values, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRows / inits.Length;
            //    if (endIndex == countRows - 1)
            //        endIndex += 1;

            //    threads[i] = new Thread(inits[i].Inititalize);
            //    threads[i].Start();
            //}

            //foreach (var thread in threads)
            //    thread.Join();
            #endregion

            #region Новая реализация потоков v2 (Parallel.For)
            Parallel.For(0, points.Count, (i) => 
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(points[i].X, j);
                }
                values[i, 0] = points[i].Y;
            });
            #endregion

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), argumentsMatrix).GetInverseAsParallel();
            Matrix second = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), valuesMatrix);
            Matrix A = Matrix.MultiplicationAsParallel(first, second);

            var count = points.Last().X / h;

            for (double x = points.First().X; x <= points.Last().X; x += h)
            {
                double y = 0;
                for (int k = 0; k <= degree; k++)
                {
                    y += A[k, 0] * Math.Pow(x, k);
                }
                resultPoints.Add(new Point(x, y));
            }

            return resultPoints;
        }

        public static List<Point> ParallelMethodOfMinimumRoots(Point[] points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Length, degree + 1]; // X
            double[,] values = new double[points.Length, 1]; // Y

            #region Старая реализация потоков
            //int countRowsArgs = arguments.GetLength(0);
            //int processorCount = Environment.ProcessorCount;

            //Initialization[] initsArgs = new Initialization[processorCount % 2 == 0 ? (processorCount / 2) : (processorCount / 2 + 1)];
            //Thread[] threadsArgs = new Thread[initsArgs.Length];

            //int startIndex = 0;
            //int endIndex = countRowsArgs / initsArgs.Length;
            //for (int i = 0; i < initsArgs.Length; i++)
            //{
            //    initsArgs[i] = new Initialization(arguments, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRowsArgs / initsArgs.Length;
            //    if (endIndex == countRowsArgs - 1)
            //        endIndex += 1;

            //    threadsArgs[i] = new Thread(initsArgs[i].InitializeArguments);
            //    threadsArgs[i].Start();
            //}

            //int countRowsValues = values.GetLength(0);
            //Initialization[] initsValues = new Initialization[processorCount / 2];
            //Thread[] threadsValues = new Thread[initsValues.Length];

            //startIndex = 0;
            //endIndex = countRowsValues / initsValues.Length;

            //for (int i = 0; i < initsValues.Length; i++)
            //{
            //    initsValues[i] = new Initialization(values, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRowsValues / initsValues.Length;
            //    if (endIndex == countRowsArgs - 1)
            //        endIndex += 1;

            //    threadsValues[i] = new Thread(initsValues[i].InitializeValues);
            //    threadsValues[i].Start();
            //}

            //foreach (var thread in threadsArgs.Concat(threadsValues))
            //{
            //    thread.Join();
            //}
            #endregion

            #region Новая реализация потоков
            //Initialization[] inits = new Initialization[Environment.ProcessorCount];
            //Thread[] threads = new Thread[Environment.ProcessorCount];

            //int countRows = points.Count;
            //int startIndex = 0;
            //int endIndex = countRows / inits.Length;
            //for (int i = 0; i < inits.Length; i++)
            //{
            //    inits[i] = new Initialization(arguments, values, points, startIndex, endIndex);
            //    startIndex = endIndex;
            //    endIndex += countRows / inits.Length;
            //    if (endIndex == countRows - 1)
            //        endIndex += 1;

            //    threads[i] = new Thread(inits[i].Inititalize);
            //    threads[i].Start();
            //}

            //foreach (var thread in threads)
            //    thread.Join();
            #endregion

            #region Новая реализация потоков v2 (Parallel.For)
            Parallel.For(0, points.Length, (i) =>
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(points[i].X, j);
                }
                values[i, 0] = points[i].Y;
            });
            #endregion

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), argumentsMatrix).GetInverseAsParallel();
            Matrix second = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), valuesMatrix);
            Matrix A = Matrix.MultiplicationAsParallel(first, second);

            for (double x = points[0].X; x <= points[points.Length - 1].X; x += h)
            {
                double y = 0;
                for (int k = 0; k <= degree; k++)
                {
                    y += A[k, 0] * Math.Pow(x, k);
                }
                resultPoints.Add(new Point(x, y));
            }

            return resultPoints;
        }

        private class Initialization
        {
            double[,] args;
            double[,] values;
            List<Point> points;
            int startIndexRow, endIndex;

            public Initialization(double[,] args, List<Point> points, int startIndexRow, int endIndex)
            {
                this.args = args;
                this.points = points;
                this.startIndexRow = startIndexRow;
                this.endIndex = endIndex;
            }

            public Initialization(double[,] args, double[,] values, List<Point> points, int startIndexRow, int endIndex)
            {
                int countRows = points.Count;

                if (values.GetLength(0) != countRows && args.GetLength(0) != countRows)
                    throw new Exception();

                if (startIndexRow < 0 || this.endIndex < 0 || startIndexRow >= countRows || this.endIndex > countRows)
                    throw new Exception();

                this.args = args;
                this.values = values;
                this.points = points;
                this.startIndexRow = startIndexRow;
                this.endIndex = endIndex;
            }

            public void Inititalize() 
            {
                for (int i = startIndexRow; i < endIndex; i++)
                {
                    for (int j = 0; j < args.GetLength(1); j++)
                    {
                        args[i, j] = Math.Pow(points[i].X, j);
                    }
                    for (int j = 0; j < values.GetLength(1); j++)
                    {
                        values[i, j] = points[i].Y;
                    }
                }
            }

            public void InitializeArguments() 
            {
                int countRows = args.GetLength(0);
                if (startIndexRow >= countRows || this.endIndex > countRows || startIndexRow < 0 || this.endIndex < 0) return;

                for (int i = startIndexRow; i < this.endIndex; i++)
                {
                    for (int j = 0; j < args.GetLength(1); j++)
                    {
                        args[i, j] = Math.Pow(points[i].X, j);
                    }
                }
            }

            public void InitializeValues()
            {
                int countRows = args.GetLength(0);
                if (startIndexRow >= countRows || this.endIndex >= countRows || startIndexRow < 0 || this.endIndex < 0) return;

                for (int i = 0; i < args.GetLength(0); i++)
                {
                    for (int j = 0; j < args.GetLength(1); j++)
                    {
                        args[i, j] = points[i].Y;
                    }
                }
            }
        }
    }
}
