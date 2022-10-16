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

            int degree = 5;
            double step = 0.001;
            List<Point> pointsList = Generator.GeneratePoints(1, 30, 1000000);
            //var pointsList = FileService.ReadFromJson<List<Point>>("data.json");

            //pointsList = pointsList.AsParallel().OrderBy(p => p.X).ToList();
            //pointsList = pointsList.OrderBy(p => p.X).ToList();
            Point[] pointsArray = pointsList.ToArray();

            //Console.WriteLine("Input data:");
            //PrintPoints(pointsArray);

            // Многопоточный метод
            stopWatch.Start();
            var resultParallel = Approximation.ParallelMethodOfMinimumRoots(pointsList, degree, step);
            stopWatch.Stop();
            Console.WriteLine("Multithreaded method time (ms): " + stopWatch.ElapsedMilliseconds);

            // Однопоточный метод
            stopWatch.Restart();
            var resultSingle = Approximation.MethodOfMinimumRoots(pointsList, degree, step);
            stopWatch.Stop();
            Console.WriteLine("Single thread method time (ms): " + stopWatch.ElapsedMilliseconds);

            // Вывод результата
            //Console.WriteLine("\nMany threads result:");
            //foreach (var p in resultParallel)
            //    Console.WriteLine($"X: {p.X}; Y: {p.Y}");

            //Console.WriteLine("\nSingle thread result:");
            //foreach (var p in resultSingle)
            //    Console.WriteLine($"X: {p.X}; Y: {p.Y}");

            Console.Read();
        }

        static void PrintPoints(IEnumerable<Point> points) 
        {
            foreach (var p in points)
                Console.WriteLine($"X: {p.X}; Y: {p.Y}");
        }
    }
}
