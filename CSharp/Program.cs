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

            int degree;
            double step;
            List<Point> pointsList;
            
            LoadDataFromJson(out degree, out step, out pointsList);
            pointsList = pointsList.OrderBy(p => p.X).ToList();

            Point[] pointsArray = pointsList.ToArray();
            HashSet<Point> pointsHashSet = new HashSet<Point>(pointsList);
            //Hashtable pointsHashtable = new Hashtable(pointsList.ToDictionary(x => x.X, y => y.Y));

            int n = 100;
            long totalTime = 0;
            for (int i = 0; i < n; i++)
            {
                stopWatch.Start();
                //Approximation.MethodOfMinimumRoots(pointsArray, degree, step); // array
                //Approximation.MethodOfMinimumRoots(pointsList, degree, step); // list
                Approximation.MethodOfMinimumRoots(pointsHashSet, degree, step); // hashset
                stopWatch.Stop();
                
                totalTime += stopWatch.ElapsedMilliseconds;
                stopWatch.Reset();
            }
            Console.WriteLine("Average time: " + totalTime / n);

            Console.Read();
        }

        static void LoadDataFromJson(out int degree, out double step, out List<Point> points)
        {
            string pathDebug = Path.GetDirectoryName(Assembly.GetExecutingAssembly().Location);
            string pathProject = Path.GetDirectoryName(Path.GetDirectoryName(Path.GetDirectoryName(pathDebug)));
            string pathJsonDataFile = Path.Combine(pathProject, "testData.json");

            JObject dataJson = JObject.Parse(File.ReadAllText(pathJsonDataFile));
            degree = dataJson.Value<int>("degree");
            step = dataJson.Value<double>("step");

            points = new List<Point>();
            JToken point = dataJson.GetValue("points").First;
            do
            {
                double x = point.Value<double>("X");
                double y = point.Value<double>("Y");
                points.Add(new Point(x, y));
                point = point.Next;
            }
            while (point != null);
        }
    }
}
