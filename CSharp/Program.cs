using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
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
            int degree;
            double step;
            List<Point> points;
            LoadDataFromJson(out degree, out step, out points);

            Approximation approximation = new Approximation(points);
            List<Point> result = approximation.MethodOfMinimumRoots(degree, step);

            foreach (Point p in result)
            {
                Console.WriteLine("X: {0}; Y: {1}", p.X, p.Y);
            }

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
