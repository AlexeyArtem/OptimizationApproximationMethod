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

            for (int i = 0; i < arguments.GetLength(0); i++)
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(points[i].X, j);
                }
            }

            for (int i = 0; i < values.GetLength(0); i++)
            {
                for (int j = 0; j < values.GetLength(1); j++)
                {
                    values[i, j] = points[i].Y;
                }
            }

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = (argumentsMatrix.GetTransporse() * argumentsMatrix).GetInverse();
            Matrix second = argumentsMatrix.GetTransporse() * valuesMatrix;
            Matrix A = first * second;

            for (double x = points[0].X; x <= points[points.Count - 1].X; x += h)
            {
                double y = 0;
                for (int i = 0; i <= degree; i++)
                {
                    y += A[i, 0] * Math.Pow(x, i);
                }
                resultPoints.Add(new Point(x, y));
            }

            return resultPoints;
        }

        public static List<Point> MethodOfMinimumRoots(Point[] points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            int sizeResult = (int)(points[points.Length - 1].X / h); // ??
            List<Point> resultPoints = new List<Point>();

            // XA=Y
            double[,] arguments = new double[points.Length, degree + 1]; // X
            double[,] values = new double[points.Length, 1]; // Y

            for (int i = 0; i < arguments.GetLength(0); i++)
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(points[i].X, j);
                }
            }

            for (int i = 0; i < values.GetLength(0); i++)
            {
                for (int j = 0; j < values.GetLength(1); j++)
                {
                    values[i, j] = points[i].Y;
                }
            }

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = (argumentsMatrix.GetTransporse() * argumentsMatrix).GetInverse();
            Matrix second = argumentsMatrix.GetTransporse() * valuesMatrix;
            Matrix A = first * second;

            int k = 0;
            for (double x = points[0].X; x <= points[points.Length - 1].X; x += h)
            {
                double y = 0;
                for (int i = 0; i <= degree; i++)
                {
                    y += A[i, 0] * Math.Pow(x, i);
                }
                resultPoints.Add(new Point(x, y));
                k++;
            }

            return resultPoints;
        }

        public static List<Point> MethodOfMinimumRoots(HashSet<Point> points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Count, degree + 1]; // X
            double[,] values = new double[points.Count, 1]; // Y

            int i = 0;
            foreach (Point p in points)
            {
                for (int j = 0; j < arguments.GetLength(1); j++)
                {
                    arguments[i, j] = Math.Pow(p.X, j);
                }
                for (int k = 0; k < values.GetLength(1); k++)
                {
                    values[i, k] = p.Y;
                }
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

        public static List<Point> ParallelMethodOfMinimumRoots(List<Point> points, int degree, double h)
        {
            if (degree <= 0) throw new Exception("Степень должна быть положительным числом");

            List<Point> resultPoints = new List<Point>();
            // XA=Y
            double[,] arguments = new double[points.Count, degree + 1]; // X
            double[,] values = new double[points.Count, 1]; // Y

            int countRowsArgs = arguments.GetLength(0);
            int processorCount = Environment.ProcessorCount;

            Initialization[] initsArgs = new Initialization[processorCount % 2 == 0 ? (processorCount / 2) : (processorCount / 2 + 1)];
            Thread[] threadsArgs = new Thread[initsArgs.Length];

            int startIndex = 0;
            int endIndex = countRowsArgs / initsArgs.Length;
            for (int i = 0; i < initsArgs.Length; i++)
            {
                initsArgs[i] = new Initialization(arguments, points, startIndex, endIndex);
                startIndex = endIndex;
                endIndex += countRowsArgs / initsArgs.Length;
                if (endIndex == countRowsArgs - 1)
                    endIndex += 1;

                threadsArgs[i] = new Thread(initsArgs[i].InitializeArguments);
                threadsArgs[i].Start();
            }

            int countRowsValues = values.GetLength(0);
            Initialization[] initsValues = new Initialization[processorCount / 2];
            Thread[] threadsValues = new Thread[initsValues.Length];

            startIndex = 0;
            endIndex = countRowsValues / initsValues.Length;

            for (int i = 0; i < initsValues.Length; i++)
            {
                initsValues[i] = new Initialization(values, points, startIndex, endIndex);
                startIndex = endIndex;
                endIndex += countRowsValues / initsValues.Length;
                if (endIndex == countRowsArgs - 1)
                    endIndex += 1;

                threadsValues[i] = new Thread(initsValues[i].InitializeValues);
                threadsValues[i].Start();
            }

            foreach (var thread in threadsArgs.Concat(threadsValues))
            {
                thread.Join();
            }

            Matrix argumentsMatrix = new Matrix(arguments);
            Matrix valuesMatrix = new Matrix(values);

            Matrix first = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), argumentsMatrix).GetTransporseAsParallel();
            Matrix second = Matrix.MultiplicationAsParallel(argumentsMatrix.GetTransporseAsParallel(), valuesMatrix);
            Matrix A = Matrix.MultiplicationAsParallel(first, second);

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

        public class Initialization
        {
            double[,] array;
            List<Point> points;
            int startIndexRow, endIndex;

            public Initialization(double[,] args, List<Point> points, int startIndexRow, int endIndex)
            {
                this.array = args;
                this.points = points;
                this.startIndexRow = startIndexRow;
                this.endIndex = endIndex;
            }

            public void InitializeArguments() 
            {
                int countRows = array.GetLength(0);
                if (startIndexRow >= countRows || this.endIndex > countRows || startIndexRow < 0 || this.endIndex < 0) return;

                for (int i = startIndexRow; i < this.endIndex; i++)
                {
                    for (int j = 0; j < array.GetLength(1); j++)
                    {
                        array[i, j] = Math.Pow(points[i].X, j);
                    }
                }
            }

            public void InitializeValues()
            {
                int countRows = array.GetLength(0);
                if (startIndexRow >= countRows || this.endIndex >= countRows || startIndexRow < 0 || this.endIndex < 0) return;

                for (int i = 0; i < array.GetLength(0); i++)
                {
                    for (int j = 0; j < array.GetLength(1); j++)
                    {
                        array[i, j] = points[i].Y;
                    }
                }
            }
        }
    }
}
