using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace PolyCurve
{
    /// <summary>
    /// Represents a point in 2D Cartesian coordinate space with X and Y coordinates.
    /// Provides basic geometric operations and utility methods for point manipulation and distance calculations.
    /// </summary>
    public class Point2D
    {
        /// <summary>
        /// Gets or sets the X-coordinate (horizontal position) of the point.
        /// Represents the position along the horizontal axis in the coordinate system.
        /// </summary>
        public double X { get; set; }

        /// <summary>
        /// Gets or sets the Y-coordinate (vertical position) of the point.
        /// Represents the position along the vertical axis in the coordinate system.
        /// </summary>
        public double Y { get; set; }

        /// <summary>
        /// Initializes a new instance of the Point2D class with the specified coordinates.
        /// Creates a point at the given X and Y position in 2D space.
        /// </summary>
        /// <param name="x">The X-coordinate (horizontal position) of the point.</param>
        /// <param name="y">The Y-coordinate (vertical position) of the point.</param>
        /// <example>
        /// <code>
        /// var point = new Point2D(3.5, -2.1);
        /// Console.WriteLine($"Point at ({point.X}, {point.Y})");
        /// </code>
        /// </example>
        public Point2D(double x, double y)
        {
            X = x;
            Y = y;
        }

        /// <summary>
        /// Calculates the Euclidean distance from this point to another Point2D using the distance formula.
        /// Returns the straight-line distance between the two points in 2D space.
        /// Formula: √[(x₂-x₁)² + (y₂-y₁)²]
        /// </summary>
        /// <param name="other">The target Point2D to calculate the distance to. Must not be null.</param>
        /// <returns>The Euclidean distance between this point and the other point as a non-negative double value.</returns>
        /// <exception cref="ArgumentNullException">Thrown when the other parameter is null.</exception>
        /// <example>
        /// <code>
        /// var point1 = new Point2D(0, 0);
        /// var point2 = new Point2D(3, 4);
        /// double distance = point1.DistanceTo(point2); // Returns 5.0
        /// Console.WriteLine($"Distance: {distance}");
        /// </code>
        /// </example>
        public double DistanceTo(Point2D other)
        {
            if (other == null) throw new ArgumentNullException(nameof(other));
            return Math.Sqrt(Math.Pow(other.X - X, 2) + Math.Pow(other.Y - Y, 2));
        }

        /// <summary>
        /// Returns a string representation of the point in standard mathematical notation.
        /// Formats the point as "(X, Y)" with the current coordinate values.
        /// </summary>
        /// <returns>A string representation of the point in the format "(X, Y)".</returns>
        /// <example>
        /// <code>
        /// var point = new Point2D(2.5, -1.8);
        /// Console.WriteLine(point.ToString()); // Output: "(2.5, -1.8)"
        /// </code>
        /// </example>
        public override string ToString() => $"({X}, {Y})";
    }
}