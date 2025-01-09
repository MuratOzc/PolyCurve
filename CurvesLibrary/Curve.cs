using MathNet.Numerics;
using MathNet.Numerics.Integration;
using MathNet.Numerics.Statistics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CurvesLibrary
{
    public class Curve
    {
        public double[]? XValues { get; }
        public double[]? YValues { get; }
        public double XMin { get; }
        public double XMax { get; }
        public int PolynomialDegree { get; }
        public double[] PolynomialCoeffs { get; }
        public Func<double, double> PolynomialFunction { get; }

        private const int DefaultPolynomialDegree = 3;
        private const int DefaultNumberOfPoints = 100;
        private const double DefaultTolerance = 1e-6;
        private const int MaxIterations = 1000;

        /// <summary>
        /// Initializes a new instance of the Curve class by fitting a polynomial to the provided data points.
        /// <para>
        /// This constructor creates a polynomial curve that best fits the given x and y values.
        /// It uses polynomial regression to determine the coefficients of the curve.
        /// </para>
        /// </summary>
        /// <param name="xValues">An array of x-coordinates of the data points.</param>
        /// <param name="yValues">An array of y-coordinates of the data points.</param>
        /// <param name="polynomialDegree">
        /// The degree of the polynomial to fit. If null, uses the default degree (usually 3).
        /// Higher degrees can provide better fits but may lead to overfitting.
        /// </param>
        /// <param name="xMin">
        /// The minimum x-value for the curve's domain. If null, uses the minimum of xValues.
        /// </param>
        /// <param name="xMax">
        /// The maximum x-value for the curve's domain. If null, uses the maximum of xValues.
        /// </param>
        /// <exception cref="ArgumentException">
        /// Thrown when xValues and yValues are null or have different lengths.
        /// </exception>
        public Curve(double[] xValues, double[] yValues, int? polynomialDegree = null, double? xMin = null, double? xMax = null)
        {
            if (xValues == null || yValues == null || xValues.Length != yValues.Length)
                throw new ArgumentException("XValues and YValues must be non-null and of equal length.");

            XValues = xValues;
            YValues = yValues;
            XMin = xMin ?? xValues.Min();
            XMax = xMax ?? xValues.Max();
            PolynomialDegree = polynomialDegree ?? DefaultPolynomialDegree;

            PolynomialCoeffs = Fit.Polynomial(xValues, yValues, PolynomialDegree);
            PolynomialFunction = CreatePolynomialFunction(PolynomialCoeffs);
        }

        /// <summary>
        /// Initializes a new instance of the Curve class using provided polynomial coefficients.
        /// <para>
        /// This constructor creates a curve directly from the given polynomial coefficients,
        /// without fitting to data points. The coefficients should be provided in ascending
        /// order of the polynomial terms (constant term first, then x, x^2, etc.).
        /// </para>
        /// </summary>
        /// <param name="polynomialCoeffs">
        /// An array of coefficients for the polynomial, in ascending order of degree.
        /// </param>
        /// <param name="xMin">The minimum x-value for the curve's domain.</param>
        /// <param name="xMax">The maximum x-value for the curve's domain.</param>
        /// <exception cref="ArgumentNullException">
        /// Thrown when polynomialCoeffs is null.
        /// </exception>
        /// <exception cref="ArgumentException">
        /// Thrown when xMin is greater than or equal to xMax.
        /// </exception>
        public Curve(double[] polynomialCoeffs, double xMin, double xMax)
        {
            PolynomialCoeffs = polynomialCoeffs ?? throw new ArgumentNullException(nameof(polynomialCoeffs));
            PolynomialDegree = polynomialCoeffs.Length - 1;
            PolynomialFunction = CreatePolynomialFunction(PolynomialCoeffs);

            XMin = xMin;
            XMax = xMax;

            if (XMin >= XMax)
                throw new ArgumentException("XMin must be less than XMax");
        }

        private static Func<double, double> CreatePolynomialFunction(double[] coeffs) =>
            x => coeffs.Select((c, i) => c * Math.Pow(x, i)).Sum();
        /// <summary>
        /// Interpolates the y-value for a given x-value using the polynomial function.
        /// <para>
        /// This method uses the polynomial coefficients to calculate the y-value
        /// corresponding to the provided x-value within the curve's domain.
        /// </para>
        /// </summary>
        /// <param name="x">The x-value at which to interpolate the curve.</param>
        /// <returns>The interpolated y-value at the given x-value.</returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the provided x-value is outside the curve's domain (XMin to XMax).
        /// </exception>
        public double Interpolate(double x)
        {
            if (x < XMin || x > XMax)
                throw new ArgumentOutOfRangeException(nameof(x), "x must be within the curve's domain.");
            return PolynomialFunction(x);
        }

        /// <summary>
        /// Finds the maximum point on the curve within its domain.
        /// <para>
        /// This method uses numerical optimization to locate the highest point
        /// on the curve. It employs a gradient ascent approach with adaptive step size.
        /// </para>
        /// </summary>
        /// <param name="tolerance">
        /// The tolerance for convergence. The search stops when the change in x
        /// is less than this value. Default is 1e-6.
        /// </param>
        /// <returns>
        /// A tuple containing the x and y coordinates of the maximum point.
        /// </returns>
        public (double X, double Y) FindMaximum(double tolerance = DefaultTolerance)
        {
            Func<double, double> derivative = x =>
            {
                double h = Math.Max(Math.Abs(x) * 1e-8, 1e-8);
                return (PolynomialFunction(x + h) - PolynomialFunction(x - h)) / (2 * h);
            };

            double x = (XMin + XMax) / 2;
            double step = (XMax - XMin) / 10;

            for (int i = 0; i < MaxIterations; i++)
            {
                double dx = derivative(x);
                if (Math.Abs(dx) < tolerance)
                    break;

                x += step * Math.Sign(dx);
                x = Math.Clamp(x, XMin, XMax);
                step *= 0.9;
            }

            // If the maximum is at the edge of the domain, return the edge value
            double yMax = PolynomialFunction(x);
            double yMin = PolynomialFunction(XMin);
            double yMax2 = PolynomialFunction(XMax);

            if (yMin > yMax && yMin > yMax2)
                return (XMin, yMin);
            if (yMax2 > yMax && yMax2 > yMin)
                return (XMax, yMax2);

            return (x, yMax);
        }


        /// <summary>
        /// Finds the minimum point on the curve within its domain.
        /// <para>
        /// This method inverts the curve and uses the FindMaximum method to locate
        /// the lowest point on the original curve.
        /// </para>
        /// </summary>
        /// <param name="tolerance">
        /// The tolerance for convergence. The search stops when the change in x
        /// is less than this value. Default is 1e-6.
        /// </param>
        /// <returns>
        /// A tuple containing the x and y coordinates of the minimum point.
        /// </returns>
        public (double X, double Y) FindMinimum(double tolerance = DefaultTolerance)
        {
            // Multiply the function by -1 to find the minimum
            Func<double, double> invertedFunction = x => -PolynomialFunction(x);
            var invertedCurve = new Curve(PolynomialCoeffs.Select(c => -c).ToArray(), XMin, XMax);
            var (x, y) = invertedCurve.FindMaximum(tolerance);
            return (x, -y);
        }

        /// <summary>
        /// Finds the roots (x-intercepts) of the curve within its domain.
        /// <para>
        /// This method uses Brent's method to find roots in intervals where the function
        /// changes sign. It also checks for roots at local extrema to catch double roots.
        /// </para>
        /// </summary>
        /// <param name="tolerance">
        /// The tolerance for root-finding. Roots are considered found when the function
        /// value is within this tolerance of zero. Default is 1e-6.
        /// </param>
        /// <returns>
        /// A list of x-values where the curve intersects the x-axis.
        /// </returns>
        public List<double> FindRoots(double tolerance = DefaultTolerance)
        {
            List<double> roots = new List<double>();

            // Use Brent's method to find roots in intervals where the function changes sign
            for (double x = XMin; x < XMax; x += (XMax - XMin) / DefaultNumberOfPoints)
            {
                double nextX = Math.Min(x + (XMax - XMin) / DefaultNumberOfPoints, XMax);
                double y1 = PolynomialFunction(x);
                double y2 = PolynomialFunction(nextX);

                if (Math.Sign(y1) != Math.Sign(y2))
                {
                    try
                    {
                        double root = MathNet.Numerics.RootFinding.Brent.FindRoot(
                    PolynomialFunction,
                    x,
                    nextX,
                    tolerance,
                    MaxIterations
                );

                        if (root >= XMin && root <= XMax)
                        {
                            roots.Add(root);
                        }
                    }
                    catch (NonConvergenceException)
                    {
                        // Root finding didn't converge, skip this interval
                    }
                }
            }

            // Find local extrema as potential double roots
            var extrema = FindDerivativeRoots(tolerance);
            foreach (var x in extrema)
            {
                double y = PolynomialFunction(x);
                if (Math.Abs(y) <= tolerance)
                {
                    roots.Add(x);
                }
            }

            // Remove duplicate roots
            roots = roots.Distinct().ToList();

            return roots;
        }

        private List<double> FindDerivativeRoots(double tolerance = DefaultTolerance)
        {
            var derivative = DerivativeCurve();
            List<double> derivativeRoots = new List<double>();

            for (double x = XMin; x < XMax; x += (XMax - XMin) / DefaultNumberOfPoints)
            {
                double nextX = Math.Min(x + (XMax - XMin) / DefaultNumberOfPoints, XMax);
                double y1 = derivative.PolynomialFunction(x);
                double y2 = derivative.PolynomialFunction(nextX);

                if (Math.Sign(y1) != Math.Sign(y2))
                {
                    try
                    {
                        double root = MathNet.Numerics.RootFinding.Brent.FindRoot(
                    derivative.PolynomialFunction,
                    x,
                    nextX,
                    tolerance,
                    MaxIterations
                );

                        if (root >= XMin && root <= XMax)
                        {
                            derivativeRoots.Add(root);
                        }
                    }
                    catch (NonConvergenceException)
                    {
                        // Root finding didn't converge, skip this interval
                    }
                }
            }

            return derivativeRoots;
        }

        /// <summary>
        /// Finds the x-values corresponding to a given y-value on the curve.
        /// <para>
        /// This method creates a new curve by subtracting the given y-value from
        /// the original curve, then finds the roots of this new curve.
        /// </para>
        /// </summary>
        /// <param name="y">The y-value to find corresponding x-values for.</param>
        /// <param name="tolerance">
        /// The tolerance for root-finding. Default is 1e-6.
        /// </param>
        /// <returns>
        /// A list of x-values where the curve has the given y-value.
        /// </returns>
        public List<double> FindXForY(double y, double tolerance = DefaultTolerance)
        {
            var newCurve = this - y;
            return newCurve.FindRoots(tolerance);
        }

        /// <summary>
        /// Calculates the definite integral of the curve between two x-values.
        /// <para>
        /// This method uses numerical integration (Newton-Cotes Trapezium Rule)
        /// to compute the area under the curve between the specified x-values.
        /// </para>
        /// </summary>
        /// <param name="x1">The lower bound of integration.</param>
        /// <param name="x2">The upper bound of integration.</param>
        /// <returns>
        /// The definite integral value, representing the signed area under the curve.
        /// </returns>
        public double CalculateIntegral(double x1, double x2)
        {
            return NewtonCotesTrapeziumRule.IntegrateComposite(PolynomialFunction, x1, x2, 1000);
        }

        /// <summary>
        /// Calculates the area under the curve (or above the x-axis) between two x-values.
        /// <para>
        /// This method computes the absolute area, meaning it takes the absolute value
        /// of the function before integration. This gives the total area between the
        /// curve and the x-axis, regardless of whether the curve is above or below the axis.
        /// </para>
        /// </summary>
        /// <param name="x1">The lower bound of the area calculation.</param>
        /// <param name="x2">The upper bound of the area calculation.</param>
        /// <returns>
        /// The area under the curve between x1 and x2, always non-negative.
        /// </returns>
        public double CalculateArea(double x1, double x2)
        {
            // The integral of the absolute value of the polynomial function
            Func<double, double> absPolynomialFunction = x => Math.Abs(PolynomialFunction(x));
            return NewtonCotesTrapeziumRule.IntegrateComposite(absPolynomialFunction, x1, x2, 1000);
        }

        /// <summary>
        /// Creates a new Curve object representing the derivative of this curve.
        /// <para>
        /// The derivative curve is computed by applying the power rule of differentiation
        /// to each term of the polynomial. The degree of the resulting polynomial is
        /// one less than the original.
        /// </para>
        /// </summary>
        /// <returns>
        /// A new Curve object representing the derivative of the current curve.
        /// </returns>
        public Curve DerivativeCurve()
        {
            var derivativeCoeffs = new double[PolynomialCoeffs.Length - 1];
            for (int i = 1; i < PolynomialCoeffs.Length; i++)
            {
                derivativeCoeffs[i - 1] = i * PolynomialCoeffs[i];
            }

            return new Curve(derivativeCoeffs, XMin, XMax);
        }

        /// <summary>
        /// Creates a new Curve object representing the integral of this curve.
        /// <para>
        /// The integral curve is computed by applying the power rule of integration
        /// to each term of the polynomial. The degree of the resulting polynomial is
        /// one more than the original.
        /// </para>
        /// </summary>
        /// <param name="c">
        /// The constant of integration. This becomes the constant term in the
        /// resulting polynomial. Default is 0.
        /// </param>
        /// <returns>
        /// A new Curve object representing the integral of the current curve.
        /// </returns>
        public Curve IntegralCurve(double c = 0)
        {
            var integralCoeffs = new double[PolynomialCoeffs.Length + 1];
            integralCoeffs[0] = c;
            for (int i = 0; i < PolynomialCoeffs.Length; i++)
            {
                integralCoeffs[i + 1] = PolynomialCoeffs[i] / (i + 1);
            }

            return new Curve(integralCoeffs, XMin, XMax);
        }

        /// <summary>
        /// Finds the inflection points of the curve within its domain.
        /// <para>
        /// Inflection points are where the curve changes from being concave upwards
        /// to concave downwards, or vice versa. They are found by locating the roots
        /// of the second derivative of the curve.
        /// </para>
        /// </summary>
        /// <param name="tolerance">
        /// The tolerance for root-finding when searching for inflection points.
        /// Default is 1e-6.
        /// </param>
        /// <returns>
        /// A list of tuples containing the x and y coordinates of inflection points.
        /// </returns>
        public List<(double X, double Y)> FindInflectionPoints(double tolerance = DefaultTolerance)
        {
            var secondDerivative = DerivativeCurve().DerivativeCurve();
            return secondDerivative.FindRoots(tolerance)
                .Select(x => (X: x, Y: PolynomialFunction(x)))
                .ToList();
        }

        /// <summary>
        /// Calculates the curvature of the curve at a given x-value.
        /// <para>
        /// Curvature is a measure of how sharply a curve bends at a given point.
        /// It is computed using the first and second derivatives of the curve.
        /// </para>
        /// </summary>
        /// <param name="x">The x-value at which to calculate the curvature.</param>
        /// <returns>
        /// The curvature value at the given x-value. A higher absolute value
        /// indicates a sharper bend in the curve.
        /// </returns>
        public double CalculateCurvature(double x)
        {
            var d1 = DerivativeCurve();
            var d2 = d1.DerivativeCurve();

            double y1 = d1.PolynomialFunction(x);
            double y2 = d2.PolynomialFunction(x);

            return Math.Abs(y2) / Math.Pow(1 + y1 * y1, 1.5);
        }

        /// <summary>
        /// Finds the local extrema (minima and maxima) of the curve within its domain.
        /// <para>
        /// Local extrema are points where the curve reaches a local minimum or maximum.
        /// They are found by locating the roots of the first derivative of the curve.
        /// </para>
        /// </summary>
        /// <param name="tolerance">
        /// The tolerance for root-finding when searching for extrema. Default is 1e-6.
        /// </param>
        /// <returns>
        /// A list of tuples containing the x and y coordinates of local extrema.
        /// </returns>
        public List<(double X, double Y)> FindLocalExtrema(double tolerance = DefaultTolerance)
        {
            var derivative = DerivativeCurve();
            return derivative.FindRoots(tolerance)
                .Select(x => (X: x, Y: PolynomialFunction(x)))
                .ToList();
        }

        /// <summary>
        /// Calculates the arc length of the curve between two x-values.
        /// <para>
        /// Arc length is the distance along the curve between two points. This method
        /// uses numerical integration to approximate the arc length.
        /// </para>
        /// </summary>
        /// <param name="x1">The starting x-value.</param>
        /// <param name="x2">The ending x-value.</param>
        /// <param name="segments">
        /// The number of segments to use for approximation. Higher values give more
        /// accurate results but take longer to compute. Default is 1000.
        /// </param>
        /// <returns>
        /// The approximate arc length of the curve between x1 and x2.
        /// </returns>
        public double CalculateArcLength(double x1, double x2, int segments = 1000)
        {
            double h = (x2 - x1) / segments;
            double sum = 0;

            for (int i = 0; i < segments; i++)
            {
                double x = x1 + i * h;
                double y1 = PolynomialFunction(x);
                double y2 = PolynomialFunction(x + h);
                sum += Math.Sqrt(h * h + (y2 - y1) * (y2 - y1));
            }

            return sum;
        }

        /// <summary>
        /// Splits the curve into two separate curves at a given x-value.
        /// <para>
        /// This method creates two new Curve objects, each representing a portion
        /// of the original curve. The split point becomes the right endpoint of the
        /// first curve and the left endpoint of the second curve.
        /// </para>
        /// </summary>
        /// <param name="x">The x-value at which to split the curve.</param>
        /// <returns>
        /// A tuple containing two new Curve objects representing the left and right
        /// parts of the original curve.
        /// </returns>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when the provided x-value is outside the curve's domain.
        /// </exception>
        public (Curve Curve1, Curve Curve2) SplitCurve(double x)
        {
            if (x <= XMin || x >= XMax)
                throw new ArgumentOutOfRangeException(nameof(x), "Split point must be within the curve's domain.");

            var leftCurve = new Curve(PolynomialCoeffs, XMin, x);
            var rightCurve = new Curve(PolynomialCoeffs, x, XMax);

            return (leftCurve, rightCurve);
        }

        /// <summary>
        /// Creates a new Curve object representing the tangent line to this curve at a given x-value.
        /// <para>
        /// The tangent line is a straight line that touches the curve at a single point,
        /// having the same slope as the curve at that point.
        /// </para>
        /// </summary>
        /// <param name="x">The x-value at which to calculate the tangent line.</param>
        /// <returns>
        /// A new Curve object representing the tangent line. This will always be
        /// a linear function (polynomial of degree 1).
        /// </returns>
        public Curve TangentLine(double x)
        {
            double a = PolynomialFunction(x);
            double b = DerivativeCurve().PolynomialFunction(x);
            double c = a - b * x;

            return new Curve(new double[] { c, b }, xMin: XMin, xMax: XMax);
        }

        /// <summary>
        /// Creates a new Curve object representing the normal line to this curve at a given x-value.
        /// <para>
        /// The normal line is a straight line perpendicular to the tangent line at the given point.
        /// It intersects the curve at a right angle.
        /// </para>
        /// </summary>
        /// <param name="x">The x-value at which to calculate the normal line.</param>
        /// <returns>
        /// A new Curve object representing the normal line. This will always be
        /// a linear function (polynomial of degree 1).
        /// </returns>
        public Curve NormalLine(double x)
        {
            double a = PolynomialFunction(x);
            double b = DerivativeCurve().PolynomialFunction(x);
            double c = a + 1 / b * x;

            return new Curve(new double[] { c, -1 / b }, xMin: XMin, xMax: XMax);
        }

        /// <summary>
        /// Creates a new Curve object representing this curve shifted horizontally.
        /// <para>
        /// Horizontal shifting moves the entire curve left or right on the x-axis.
        /// A positive shift moves the curve to the right, while a negative shift
        /// moves it to the left.
        /// </para>
        /// <para>
        /// The operators << and >> can be used as shortcuts for this method.
        /// </para>
        /// </summary>
        /// <param name="shift">The amount to shift the curve horizontally.</param>
        /// <returns>
        /// A new Curve object representing the shifted curve. The polynomial
        /// coefficients are adjusted to reflect the shift.
        /// </returns>
        public Curve ShiftHorizontal(double shift)
        {
            if (shift == 0) return this;

            double[] newCoeffs = new double[PolynomialCoeffs.Length];
            for (int n = 0; n < PolynomialCoeffs.Length; n++)
            {
                for (int k = n; k < PolynomialCoeffs.Length; k++)
                {
                    newCoeffs[n] += PolynomialCoeffs[k] * BinomialCoefficient(k, n) * Math.Pow(-shift, k - n);
                }
            }

            return new Curve(newCoeffs, XMin + shift, XMax + shift);
        }

        private static long BinomialCoefficient(int n, int k)
        {
            if (k < 0 || k > n) return 0;
            if (k == 0 || k == n) return 1;

            long result = 1;
            for (int i = 1; i <= k; i++)
            {
                result *= n - (k - i);
                result /= i;
            }

            return result;
        }


        /// <summary>
        /// Finds the intersection points between this curve and another curve.
        /// <para>
        /// This method calculates the x-coordinates where the two curves intersect
        /// by finding the roots of the difference between the two curves' functions.
        /// It then calculates the corresponding y-coordinates.
        /// </para>
        /// </summary>
        /// <param name="otherCurve">The other Curve object to find intersections with.</param>
        /// <param name="tolerance">
        /// The tolerance for root-finding when searching for intersections. 
        /// Default is 1e-6.
        /// </param>
        /// <returns>
        /// A list of tuples, each containing the x and y coordinates of an intersection point.
        /// Returns an empty list if no intersections are found.
        /// </returns>
        /// <exception cref="ArgumentNullException">
        /// Thrown when otherCurve is null.
        /// </exception>
        public List<(double X, double Y)> FindIntersections(Curve otherCurve, double tolerance = 1e-6)
        {
            if (otherCurve == null)
                throw new ArgumentNullException(nameof(otherCurve));

            // Create a new curve representing the difference between this curve and the other curve
            Curve differenceCurve = this - otherCurve;

            // Find the roots of the difference curve (these are the x-coordinates of intersections)
            List<double> intersectionXs = differenceCurve.FindRoots(tolerance);

            // Calculate the y-coordinates for each intersection point
            List<(double X, double Y)> intersectionPoints = intersectionXs
            .Select(x => (X: x, Y: this.Interpolate(x)))
            .ToList();

            return intersectionPoints;
        }



        /// <summary>
        /// Interpolates this curve with another curve based on a given ratio.
        /// <para>
        /// This method creates a new curve that represents an intermediate state
        /// between this curve and the other curve. The ratio determines how close
        /// the resulting curve is to each of the original curves.
        /// </para>
        /// </summary>
        /// <param name="otherCurve">The curve to interpolate with.</param>
        /// <param name="ratio">
        /// A value between 0 and 1 that determines the interpolation:
        /// 0 results in this curve, 1 results in the other curve, and values in
        /// between create intermediate curves.
        /// </param>
        /// <returns>
        /// A new Curve object representing the interpolated curve.
        /// </returns>
        /// <exception cref="ArgumentNullException">
        /// Thrown when otherCurve is null.
        /// </exception>
        /// <exception cref="ArgumentOutOfRangeException">
        /// Thrown when ratio is not between 0 and 1.
        /// </exception>
        public Curve InterpolateWith(Curve otherCurve, double ratio)
        {
            if (otherCurve == null)
                throw new ArgumentNullException(nameof(otherCurve));
            if (ratio < 0 || ratio > 1)
                throw new ArgumentOutOfRangeException(nameof(ratio), "Ratio must be between 0 and 1.");

            return PointBasedInterpolation(otherCurve, ratio);
        }
        private Curve PointBasedInterpolation(Curve otherCurve, double ratio)
        {
            int samplePoints = 100;
            // Create a list of 100 equally spaced points for each curve
            double[] xValues1 = Enumerable.Range(0, samplePoints)
            .Select(i => this.XMin + i * (this.XMax - this.XMin) / (samplePoints - 1))
            .ToArray();
            double[] yValues1 = xValues1.Select(x => this.PolynomialFunction(x)).ToArray();

            double[] xValues2 = Enumerable.Range(0, samplePoints)
            .Select(i => otherCurve.XMin + i * (otherCurve.XMax - otherCurve.XMin) / (samplePoints - 1))
            .ToArray();
            double[] yValues2 = xValues2.Select(x => otherCurve.PolynomialFunction(x)).ToArray();

            // Interpolate the values using the ratio
            double[] xValues = new double[samplePoints];
            double[] yValues = new double[samplePoints];

            for (int i = 0; i < samplePoints; i++)
            {
                xValues[i] = (1 - ratio) * xValues1[i] + ratio * xValues2[i];
                yValues[i] = (1 - ratio) * yValues1[i] + ratio * yValues2[i];
            }

            return new Curve(xValues, yValues);
        }

        // Curve Statistics
        public class CurveStatistics
        {
            public double Mean { get; set; }
            public double Variance { get; set; }
            public double Skewness { get; set; }
            public double Kurtosis { get; set; }
        }
        /// <summary>
        /// Calculates various statistical measures for the curve.
        /// <para>
        /// This method computes the following statistics:
        /// <list type="bullet">
        /// <item>
        ///     <term>Mean</term>
        ///     <description>The average value of the curve over its domain.</description>
        /// </item>
        /// <item>
        ///     <term>Variance</term>
        ///     <description>A measure of the spread of the curve values around the mean.</description>
        /// </item>
        /// <item>
        ///     <term>Skewness</term>
        ///     <description>A measure of the asymmetry of the curve's distribution.</description>
        /// </item>
        /// <item>
        ///     <term>Kurtosis</term>
        ///     <description>A measure of the "tailedness" of the curve's distribution.</description>
        /// </item>
        /// </list>
        /// </para>
        /// </summary>
        /// <param name="numberOfPoints">
        /// The number of sample points to use for calculations. Higher values give more
        /// accurate results but take longer to compute. Default is 100.
        /// </param>
        /// <returns>
        /// A CurveStatistics object containing the calculated mean, variance, skewness, and kurtosis.
        /// </returns>
        public CurveStatistics CalculateStatistics(int numberOfPoints = DefaultNumberOfPoints)
        {
            double[] samplePoints = GenerateSamplePoints(numberOfPoints);

            return new CurveStatistics
            {
                Mean = samplePoints.Average(),
                Variance = samplePoints.Variance(),
                Skewness = samplePoints.Skewness(),
                Kurtosis = samplePoints.Kurtosis()
            };
        }

        private double[] GenerateSamplePoints(int numberOfPoints)
        {
            double step = (XMax - XMin) / (numberOfPoints - 1);
            return Enumerable.Range(0, numberOfPoints)
                .Select(i => XMin + i * step)
                .Select(x => PolynomialFunction(x))
                .ToArray();
        }

        // Curve Fitting Quality Metrics
        public class FittingQualityMetrics
        {
            public double RSquared { get; set; }
            public double MeanSquaredError { get; set; }
            public double RootMeanSquaredError { get; set; }
            public double MeanAbsoluteError { get; set; }
        }

        /// <summary>
        /// Calculates quality metrics for the curve fitting.
        /// <para>
        /// This method computes the following metrics:
        /// <list type="bullet">
        /// <item>
        ///     <term>R-squared</term>
        ///     <description>A measure of how well the curve fits the original data points.</description>
        /// </item>
        /// <item>
        ///     <term>Mean Squared Error (MSE)</term>
        ///     <description>The average squared difference between predicted and actual y-values.</description>
        /// </item>
        /// <item>
        ///     <term>Root Mean Squared Error (RMSE)</term>
        ///     <description>The square root of the MSE, providing an error measure in the same units as the y-values.</description>
        /// </item>
        /// <item>
        ///     <term>Mean Absolute Error (MAE)</term>
        ///     <description>The average absolute difference between predicted and actual y-values.</description>
        /// </item>
        /// </list>
        /// </para>
        /// </summary>
        /// <returns>
        /// A FittingQualityMetrics object containing the calculated R-squared, MSE, RMSE, and MAE values.
        /// </returns>
        public FittingQualityMetrics CalculateFittingQualityMetrics()
        {
            if (XValues == null || YValues == null || XValues.Length != YValues.Length)
                throw new InvalidOperationException("Original data points are not available or invalid.");

            double sst = YValues.Variance() * (YValues.Length - 1);  // Total sum of squares
            double sse = 0;  // Sum of squared errors
            double sae = 0;  // Sum of absolute errors

            for (int i = 0; i < XValues.Length; i++)
            {
                double error = YValues[i] - PolynomialFunction(XValues[i]);
                sse += error * error;
                sae += Math.Abs(error);
            }

            int n = XValues.Length;
            double mse = sse / n;
            double mae = sae / n;

            return new FittingQualityMetrics
            {
                RSquared = 1 - (sse / sst),
                MeanSquaredError = mse,
                RootMeanSquaredError = Math.Sqrt(mse),
                MeanAbsoluteError = mae
            };
        }


        public static Curve operator +(Curve a, Curve b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            int maxDegree = Math.Max(a.PolynomialDegree, b.PolynomialDegree);
            double[] newCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= maxDegree; i++)
            {
                double aCoeff = i < a.PolynomialCoeffs.Length ? a.PolynomialCoeffs[i] : 0;
                double bCoeff = i < b.PolynomialCoeffs.Length ? b.PolynomialCoeffs[i] : 0;
                newCoeffs[i] = aCoeff + bCoeff;
            }

            return new Curve(newCoeffs, Math.Max(a.XMin, b.XMin), Math.Min(a.XMax, b.XMax));
        }

        public static Curve operator -(Curve a, Curve b)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (b == null) throw new ArgumentNullException(nameof(b));

            int maxDegree = Math.Max(a.PolynomialDegree, b.PolynomialDegree);
            double[] newCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= maxDegree; i++)
            {
                double aCoeff = i < a.PolynomialCoeffs.Length ? a.PolynomialCoeffs[i] : 0;
                double bCoeff = i < b.PolynomialCoeffs.Length ? b.PolynomialCoeffs[i] : 0;
                newCoeffs[i] = aCoeff - bCoeff;
            }

            return new Curve(newCoeffs, Math.Max(a.XMin, b.XMin), Math.Min(a.XMax, b.XMax));
        }

        public static Curve operator *(Curve a, double scalar)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));

            double[] newCoeffs = a.PolynomialCoeffs.Select(c => c * scalar).ToArray();
            return new Curve(newCoeffs, a.XMin, a.XMax);
        }

        public static Curve operator *(double scalar, Curve a) => a * scalar;

        public static Curve operator /(Curve a, double scalar)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));
            if (scalar == 0) throw new DivideByZeroException("Cannot divide a curve by zero.");

            double[] newCoeffs = a.PolynomialCoeffs.Select(c => c / scalar).ToArray();
            return new Curve(newCoeffs, a.XMin, a.XMax);
        }

        // New operator overloads for scalar addition and subtraction
        public static Curve operator +(Curve a, double scalar)
        {
            if (a == null) throw new ArgumentNullException(nameof(a));

            double[] newCoeffs = new double[a.PolynomialCoeffs.Length];
            Array.Copy(a.PolynomialCoeffs, newCoeffs, a.PolynomialCoeffs.Length);
            newCoeffs[0] += scalar; // Add the scalar to the constant term

            return new Curve(newCoeffs, a.XMin, a.XMax);
        }

        public static Curve operator +(double scalar, Curve a) => a + scalar;

        public static Curve operator -(Curve a, double scalar) => a + (-scalar);

        public static Curve operator -(double scalar, Curve a) => a + (-scalar);

        public static Curve operator >>(Curve curve, double shift) => curve.ShiftHorizontal(shift);
        public static Curve operator <<(Curve curve, double shift) => curve.ShiftHorizontal(-shift);


        public override string ToString()
        {
            return string.Join(" + ", PolynomialCoeffs.Select((c, i) => ($"{c}x^{i}").Replace(",", "."))
                .Where(s => !string.IsNullOrEmpty(s)));
        }
    }
}
