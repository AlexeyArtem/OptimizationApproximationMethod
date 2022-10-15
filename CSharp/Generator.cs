using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;

namespace CSharp
{
    static class Generator
    {
        public static List<Point> GeneratePoints(int minValue, int maxValue, int count)
        {
            List<Point> points = new List<Point>(count);
            Random random = new Random();

            for (int i = 0; i < count; i++)
            {
                double x = random.Next(minValue, maxValue) + Math.Round(random.NextDouble(), 2);
                double y = random.Next(minValue, maxValue) + Math.Round(random.NextDouble(), 2);
                points.Add(new Point(x, y));
            }

            return points;
        }
    }
}
