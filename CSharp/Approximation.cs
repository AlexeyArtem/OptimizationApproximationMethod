using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
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
    }
}
