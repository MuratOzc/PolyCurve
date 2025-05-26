using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PolyCurve
{

    /// <summary>
    /// Represents a point in 2D space with X and Y coordinates.
    /// </summary>
    public class Point2D
    {
        public double X { get; set; }
        public double Y { get; set; }

        public Point2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// Calculates the distance to another Point2D.
        /// </summary>
        /// <param name="other">The other Point2D to calculate the distance to.</param>
        /// <returns>The distance between the two points.</returns>
        public double DistanceTo(Point2D other)
        {
            if (other == null) throw new ArgumentNullException(nameof(other));
            return Math.Sqrt(Math.Pow(other.X - X, 2) + Math.Pow(other.Y - Y, 2));
        }

        public override string ToString() => $"({X}, {Y})";
    }
}
