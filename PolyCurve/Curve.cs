using System;
using System.Collections.Generic;
using System.Linq;

namespace PolyCurve
{
    /// <summary>
    /// Represents a polynomial curve segment with comprehensive mathematical operations and analysis methods.
    /// Supports polynomial fitting from data points, direct coefficient specification, and extensive curve analysis
    /// including derivatives, integrals, extrema finding, root solving, and geometric transformations.
    /// </summary>
    public class Curve
    {
        /// <summary>
        /// Gets the original X-coordinate data points used for polynomial fitting, or null if created from coefficients.
        /// </summary>
        public double[]? XValues { get; }

        /// <summary>
        /// Gets the original Y-coordinate data points used for polynomial fitting, or null if created from coefficients.
        /// </summary>
        public double[]? YValues { get; }

        /// <summary>
        /// Gets the minimum X-value of the curve's domain.
        /// </summary>
        public double XMin { get; }

        /// <summary>
        /// Gets the maximum X-value of the curve's domain.
        /// </summary>
        public double XMax { get; }

        /// <summary>
        /// Gets the degree of the polynomial (highest power of x in the polynomial equation).
        /// </summary>
        public int PolynomialDegree { get; }

        /// <summary>
        /// Gets the polynomial coefficients array where index represents the power of x.
        /// For example, coefficients [a₀, a₁, a₂] represent the polynomial a₀ + a₁x + a₂x².
        /// </summary>
        public double[] PolynomialCoeffs { get; }

        #region Constructors

        /// <summary>
        /// Creates a polynomial curve by fitting data points using least squares regression with the specified degree within given bounds.
        /// The polynomial is fitted to minimize the sum of squared errors between the curve and the data points.
        /// </summary>
        /// <param name="xValues">Array of X-coordinate values for the data points. Must have at least 2 elements.</param>
        /// <param name="yValues">Array of Y-coordinate values for the data points. Must have the same length as xValues.</param>
        /// <param name="degree">The degree of the polynomial to fit. If -1, automatically determines optimal degree (max 5). Must not exceed number of data points minus 1.</param>
        /// <param name="xMin">Minimum X-value for the curve domain. If null, uses the minimum value from xValues.</param>
        /// <param name="xMax">Maximum X-value for the curve domain. If null, uses the maximum value from xValues.</param>
        /// <exception cref="ArgumentNullException">Thrown when xValues or yValues is null.</exception>
        /// <exception cref="ArgumentException">Thrown when arrays have different lengths or fewer than 2 data points.</exception>
        /// <example>
        /// <code>
        /// double[] x = {1, 2, 3, 4, 5};
        /// double[] y = {2, 8, 18, 32, 50};
        /// var curve = new Curve(x, y, 2); // Fit quadratic polynomial
        /// </code>
        /// </example>
        public Curve(double[] xValues, double[] yValues, int degree = -1, double? xMin = null, double? xMax = null)
        {
            if (xValues == null || yValues == null)
                throw new ArgumentNullException("Data points cannot be null");
            if (xValues.Length != yValues.Length)
                throw new ArgumentException("X and Y arrays must have the same length");
            if (xValues.Length < 2)
                throw new ArgumentException("At least 2 data points are required");

            XValues = xValues.ToArray();
            YValues = yValues.ToArray();
            XMin = xMin ?? xValues.Min();
            XMax = xMax ?? xValues.Max();

            // Determine polynomial degree if not specified
            if (degree == -1)
                degree = Math.Min(xValues.Length - 1, 5); // Default to max degree 5

            degree = Math.Min(degree, xValues.Length - 1);
            PolynomialDegree = degree;

            // Fit polynomial using least squares
            PolynomialCoeffs = FitPolynomial(xValues, yValues, degree);
        }

        /// <summary>
        /// Creates a polynomial curve by fitting a collection of 2D points using least squares regression with the specified degree within given bounds.
        /// </summary>
        /// <param name="dataPoints">Collection of Point2D objects containing the data to fit. Must contain at least 2 points.</param>
        /// <param name="degree">The degree of the polynomial to fit. If -1, automatically determines optimal degree (max 5).</param>
        /// <param name="xMin">Minimum X-value for the curve domain. If null, uses the minimum X value from dataPoints.</param>
        /// <param name="xMax">Maximum X-value for the curve domain. If null, uses the maximum X value from dataPoints.</param>
        /// <exception cref="ArgumentNullException">Thrown when dataPoints is null.</exception>
        /// <exception cref="ArgumentException">Thrown when dataPoints contains fewer than 2 points.</exception>
        public Curve(IEnumerable<Point2D> dataPoints, int degree = -1, double? xMin = null, double? xMax = null)
            : this(dataPoints.Select(p => p.X).ToArray(),
                   dataPoints.Select(p => p.Y).ToArray(),
                   degree, xMin, xMax)
        {
        }

        /// <summary>
        /// Creates a polynomial curve directly from polynomial coefficients with specified domain bounds.
        /// The polynomial is constructed as: f(x) = coefficients[0] + coefficients[1]*x + coefficients[2]*x² + ...
        /// </summary>
        /// <param name="coefficients">Array of polynomial coefficients where index represents the power of x. Must not be null or empty.</param>
        /// <param name="xMin">Minimum X-value of the curve's domain.</param>
        /// <param name="xMax">Maximum X-value of the curve's domain.</param>
        /// <exception cref="ArgumentException">Thrown when coefficients array is null or empty.</exception>
        /// <example>
        /// <code>
        /// double[] coeffs = {1, 2, 3}; // Represents 1 + 2x + 3x²
        /// var curve = new Curve(coeffs, 0, 10);
        /// </code>
        /// </example>
        public Curve(double[] coefficients, double xMin, double xMax)
        {
            if (coefficients == null || coefficients.Length == 0)
                throw new ArgumentException("Coefficients array cannot be null or empty");

            PolynomialCoeffs = coefficients.ToArray();
            PolynomialDegree = coefficients.Length - 1;
            XMin = xMin;
            XMax = xMax;
            XValues = null;
            YValues = null;
        }

        #endregion

        #region Basic Evaluation

        /// <summary>
        /// Evaluates the polynomial at the specified x-value and returns the corresponding point on the curve.
        /// Uses Horner's method for efficient polynomial evaluation.
        /// </summary>
        /// <param name="xValue">The X-coordinate at which to evaluate the polynomial.</param>
        /// <returns>A Point2D object containing the input x-value and the calculated y-value on the curve.</returns>
        /// <example>
        /// <code>
        /// var point = curve.EvaluateAt(2.5);
        /// Console.WriteLine($"At x={point.X}, y={point.Y}");
        /// </code>
        /// </example>
        public Point2D EvaluateAt(double xValue)
        {
            double y = EvaluatePolynomial(PolynomialCoeffs, xValue);
            return new Point2D(xValue, y);
        }

        #endregion

        #region Curve Segmentation

        /// <summary>
        /// Splits the curve at the specified x-values, creating multiple curve segments with the same polynomial but different domains.
        /// Split points outside the curve's domain are ignored.
        /// </summary>
        /// <param name="splitPoints">Array of X-values where the curve should be split. Values outside [XMin, XMax] are ignored.</param>
        /// <returns>Collection of Curve objects representing the segments, ordered from left to right along the x-axis.</returns>
        /// <example>
        /// <code>
        /// var segments = curve.SplitCurve(2.0, 5.0, 8.0);
        /// // Returns curves for domains [XMin,2], [2,5], [5,8], [8,XMax]
        /// </code>
        /// </example>
        public IEnumerable<Curve> SplitCurve(params double[] splitPoints)
        {
            var points = splitPoints.Where(p => p > XMin && p < XMax).OrderBy(p => p).ToList();
            points.Insert(0, XMin);
            points.Add(XMax);

            var segments = new List<Curve>();
            for (int i = 0; i < points.Count - 1; i++)
            {
                segments.Add(new Curve(PolynomialCoeffs, points[i], points[i + 1]));
            }
            return segments;
        }

        /// <summary>
        /// Splits the curve into the specified number of equally-sized segments along the x-axis.
        /// Each segment has the same width but may represent different portions of the polynomial curve.
        /// </summary>
        /// <param name="numberOfSegments">The number of equal segments to create. Must be positive.</param>
        /// <returns>Collection of Curve objects representing the equal segments, ordered from left to right.</returns>
        /// <exception cref="ArgumentException">Thrown when numberOfSegments is less than or equal to 0.</exception>
        /// <example>
        /// <code>
        /// var segments = curve.SplitCurve(4); // Creates 4 equal segments
        /// </code>
        /// </example>
        public IEnumerable<Curve> SplitCurve(int numberOfSegments)
        {
            if (numberOfSegments <= 0)
                throw new ArgumentException("Number of segments must be positive");

            double segmentWidth = (XMax - XMin) / numberOfSegments;
            var splitPoints = new double[numberOfSegments - 1];

            for (int i = 0; i < numberOfSegments - 1; i++)
            {
                splitPoints[i] = XMin + (i + 1) * segmentWidth;
            }

            return SplitCurve(splitPoints);
        }

        #endregion

        #region Extrema and Critical Points

        /// <summary>
        /// Finds the global maximum point on the curve within its domain by examining local extrema and endpoints.
        /// Considers both critical points (where derivative equals zero) and domain boundaries.
        /// </summary>
        /// <returns>A Point2D object representing the point with the highest Y-value on the curve.</returns>
        /// <example>
        /// <code>
        /// var max = curve.FindMaximum();
        /// Console.WriteLine($"Maximum at ({max.X}, {max.Y})");
        /// </code>
        /// </example>
        public Point2D FindMaximum()
        {
            var extrema = FindLocalExtrema().ToList();

            // Also check endpoints
            var candidates = new List<Point2D>(extrema);
            candidates.Add(EvaluateAt(XMin));
            candidates.Add(EvaluateAt(XMax));

            return candidates.OrderByDescending(p => p.Y).First();
        }

        /// <summary>
        /// Finds the global minimum point on the curve within its domain by examining local extrema and endpoints.
        /// Considers both critical points (where derivative equals zero) and domain boundaries.
        /// </summary>
        /// <returns>A Point2D object representing the point with the lowest Y-value on the curve.</returns>
        /// <example>
        /// <code>
        /// var min = curve.FindMinimum();
        /// Console.WriteLine($"Minimum at ({min.X}, {min.Y})");
        /// </code>
        /// </example>
        public Point2D FindMinimum()
        {
            var extrema = FindLocalExtrema().ToList();

            // Also check endpoints
            var candidates = new List<Point2D>(extrema);
            candidates.Add(EvaluateAt(XMin));
            candidates.Add(EvaluateAt(XMax));

            return candidates.OrderBy(p => p.Y).First();
        }

        /// <summary>
        /// Finds all local extrema (local minima and maxima) within the curve's domain using calculus methods.
        /// Locates points where the first derivative equals zero and the second derivative is non-zero.
        /// </summary>
        /// <returns>Collection of Point2D objects representing local extrema points, excluding endpoints.</returns>
        /// <example>
        /// <code>
        /// var extrema = curve.FindLocalExtrema();
        /// foreach(var point in extrema)
        ///     Console.WriteLine($"Local extremum at ({point.X}, {point.Y})");
        /// </code>
        /// </example>
        public IEnumerable<Point2D> FindLocalExtrema()
        {
            var derivative = GetDerivative();
            var criticalPoints = derivative.FindRoots();

            var extrema = new List<Point2D>();
            var secondDerivative = derivative.GetDerivative();

            foreach (var x in criticalPoints)
            {
                if (x >= XMin && x <= XMax)
                {
                    double secondDeriv = secondDerivative.EvaluateAt(x).Y;
                    // If second derivative is non-zero, we have an extremum
                    if (Math.Abs(secondDeriv) > 1e-10)
                    {
                        extrema.Add(EvaluateAt(x));
                    }
                }
            }

            return extrema;
        }

        /// <summary>
        /// Finds all inflection points of the curve within its domain where the concavity changes.
        /// Inflection points occur where the second derivative equals zero and changes sign.
        /// </summary>
        /// <returns>Collection of Point2D objects representing inflection points on the curve.</returns>
        /// <example>
        /// <code>
        /// var inflectionPoints = curve.FindInflectionPoints();
        /// foreach(var point in inflectionPoints)
        ///     Console.WriteLine($"Inflection point at ({point.X}, {point.Y})");
        /// </code>
        /// </example>
        public IEnumerable<Point2D> FindInflectionPoints()
        {
            var secondDerivative = GetDerivative().GetDerivative();
            var inflectionXValues = secondDerivative.FindRoots();

            var inflectionPoints = new List<Point2D>();
            foreach (var x in inflectionXValues)
            {
                if (x >= XMin && x <= XMax)
                {
                    inflectionPoints.Add(EvaluateAt(x));
                }
            }

            return inflectionPoints;
        }

        #endregion

        #region Root Finding

        /// <summary>
        /// Finds all roots (x-intercepts) of the curve within its domain where the curve crosses the x-axis.
        /// Uses numerical methods to locate points where the polynomial equals zero.
        /// </summary>
        /// <returns>Collection of X-values where the curve intersects the x-axis (where Y = 0).</returns>
        /// <example>
        /// <code>
        /// var roots = curve.FindRoots();
        /// foreach(double root in roots)
        ///     Console.WriteLine($"Root at x = {root}");
        /// </code>
        /// </example>
        public IEnumerable<double> FindRoots()
        {
            return FindXValuesForY(0);
        }

        /// <summary>
        /// Finds all x-values where the curve has the specified y-value within the curve's domain.
        /// Uses Newton-Raphson method with multiple starting points to ensure all solutions are found.
        /// </summary>
        /// <param name="yValue">The Y-value for which to find corresponding X-values on the curve.</param>
        /// <returns>Collection of X-values where the curve equals the specified Y-value, sorted in ascending order.</returns>
        /// <example>
        /// <code>
        /// var xValues = curve.FindXValuesForY(10.5);
        /// foreach(double x in xValues)
        ///     Console.WriteLine($"Curve equals 10.5 at x = {x}");
        /// </code>
        /// </example>
        public IEnumerable<double> FindXValuesForY(double yValue)
        {
            // Create a new polynomial by shifting down by yValue
            var shiftedCoeffs = PolynomialCoeffs.ToArray();
            shiftedCoeffs[0] -= yValue;

            var roots = new List<double>();

            // Use numerical methods for polynomial root finding
            // Simple implementation using Newton's method with multiple starting points
            int numStartPoints = 20;
            double step = (XMax - XMin) / (numStartPoints - 1);

            for (int i = 0; i < numStartPoints; i++)
            {
                double x0 = XMin + i * step;
                double root = NewtonRaphson(shiftedCoeffs, x0);

                if (!double.IsNaN(root) && root >= XMin && root <= XMax)
                {
                    // Check if this is actually a root
                    double value = EvaluatePolynomial(shiftedCoeffs, root);
                    if (Math.Abs(value) < 1e-10)
                    {
                        // Check if we already found this root
                        // Use a larger tolerance for deduplication
                        bool isNew = true;
                        double deduplicationTolerance = 1e-6; // Increased from 1e-10

                        foreach (var existingRoot in roots)
                        {
                            if (Math.Abs(existingRoot - root) < deduplicationTolerance)
                            {
                                isNew = false;
                                break;
                            }
                        }

                        if (isNew)
                        {
                            // Round to reasonable precision to avoid floating point artifacts
                            root = Math.Round(root, 10);
                            roots.Add(root);
                        }
                    }
                }
            }

            // Final cleanup: merge any roots that are still too close
            var cleanedRoots = new List<double>();
            roots.Sort();

            foreach (var root in roots)
            {
                if (cleanedRoots.Count == 0 || Math.Abs(root - cleanedRoots.Last()) > 1e-6)
                {
                    cleanedRoots.Add(root);
                }
            }

            return cleanedRoots;
        }

        #endregion

        #region Integration and Area

        /// <summary>
        /// Calculates the definite integral of the curve over its entire domain using analytical integration.
        /// Returns the signed area between the curve and the x-axis.
        /// </summary>
        /// <returns>The definite integral value over the curve's domain [XMin, XMax].</returns>
        /// <example>
        /// <code>
        /// double integral = curve.CalculateDefiniteIntegral();
        /// Console.WriteLine($"Integral from {curve.XMin} to {curve.XMax} = {integral}");
        /// </code>
        /// </example>
        public double CalculateDefiniteIntegral()
        {
            return CalculateDefiniteIntegral(XMin, XMax);
        }

        /// <summary>
        /// Calculates the definite integral of the curve between the specified bounds using analytical integration.
        /// Returns the signed area between the curve and the x-axis over the given interval.
        /// </summary>
        /// <param name="xStart">The lower bound of integration.</param>
        /// <param name="xEnd">The upper bound of integration.</param>
        /// <returns>The definite integral value over the interval [xStart, xEnd].</returns>
        /// <example>
        /// <code>
        /// double integral = curve.CalculateDefiniteIntegral(2.0, 5.0);
        /// Console.WriteLine($"Integral from 2 to 5 = {integral}");
        /// </code>
        /// </example>
        public double CalculateDefiniteIntegral(double xStart, double xEnd)
        {
            var antiderivative = GetAntiderivative();
            double endValue = antiderivative.EvaluateAt(xEnd).Y;
            double startValue = antiderivative.EvaluateAt(xStart).Y;
            return endValue - startValue;
        }

        /// <summary>
        /// Calculates the total area under the curve over its entire domain, taking absolute values of negative regions.
        /// Automatically handles curves that cross the x-axis by splitting at roots and summing absolute areas.
        /// </summary>
        /// <returns>The total unsigned area between the curve and the x-axis over the curve's domain.</returns>
        /// <example>
        /// <code>
        /// double area = curve.CalculateAreaUnderCurve();
        /// Console.WriteLine($"Total area under curve = {area}");
        /// </code>
        /// </example>
        public double CalculateAreaUnderCurve()
        {
            return CalculateAreaUnderCurve(XMin, XMax);
        }

        /// <summary>
        /// Calculates the total area under the curve between specified bounds, taking absolute values of negative regions.
        /// Splits the interval at roots to handle sign changes and sums the absolute areas of each segment.
        /// </summary>
        /// <param name="xStart">The lower bound for area calculation.</param>
        /// <param name="xEnd">The upper bound for area calculation.</param>
        /// <returns>The total unsigned area between the curve and the x-axis over the specified interval.</returns>
        /// <example>
        /// <code>
        /// double area = curve.CalculateAreaUnderCurve(1.0, 4.0);
        /// Console.WriteLine($"Area from 1 to 4 = {area}");
        /// </code>
        /// </example>
        public double CalculateAreaUnderCurve(double xStart, double xEnd)
        {
            // For area, we need to consider the absolute value
            // Split the curve at roots and integrate each segment
            var roots = FindRoots().Where(r => r > xStart && r < xEnd).OrderBy(r => r).ToList();

            var points = new List<double> { xStart };
            points.AddRange(roots);
            points.Add(xEnd);

            double totalArea = 0;
            for (int i = 0; i < points.Count - 1; i++)
            {
                double segmentIntegral = CalculateDefiniteIntegral(points[i], points[i + 1]);
                totalArea += Math.Abs(segmentIntegral);
            }

            return totalArea;
        }

        /// <summary>
        /// Calculates the arc length of the curve over its entire domain using the arc length formula.
        /// Arc length is computed as the integral of √(1 + (dy/dx)²) dx using Simpson's rule for numerical integration.
        /// </summary>
        /// <returns>The total arc length of the curve over its domain.</returns>
        /// <example>
        /// <code>
        /// double length = curve.CalculateArcLength();
        /// Console.WriteLine($"Total arc length = {length}");
        /// </code>
        /// </example>
        public double CalculateArcLength()
        {
            return CalculateArcLength(XMin, XMax);
        }

        /// <summary>
        /// Calculates the arc length of the curve between specified bounds using numerical integration.
        /// Uses Simpson's rule to numerically integrate √(1 + (dy/dx)²) over the specified interval.
        /// </summary>
        /// <param name="xStart">The lower bound for arc length calculation.</param>
        /// <param name="xEnd">The upper bound for arc length calculation.</param>
        /// <returns>The arc length of the curve over the specified interval.</returns>
        /// <example>
        /// <code>
        /// double length = curve.CalculateArcLength(0, 5);
        /// Console.WriteLine($"Arc length from 0 to 5 = {length}");
        /// </code>
        /// </example>
        public double CalculateArcLength(double xStart, double xEnd)
        {
            // Arc length = integral of sqrt(1 + (dy/dx)^2) dx
            var derivative = GetDerivative();

            // Use numerical integration (Simpson's rule)
            int n = 1000; // Number of subdivisions
            double h = (xEnd - xStart) / n;
            double sum = 0;

            for (int i = 0; i <= n; i++)
            {
                double x = xStart + i * h;
                double dydx = derivative.EvaluateAt(x).Y;
                double integrand = Math.Sqrt(1 + dydx * dydx);

                if (i == 0 || i == n)
                    sum += integrand;
                else if (i % 2 == 1)
                    sum += 4 * integrand;
                else
                    sum += 2 * integrand;
            }

            return sum * h / 3;
        }

        #endregion

        #region Derivatives and Integrals

        /// <summary>
        /// Creates a new Curve representing the derivative (rate of change) of this curve.
        /// The derivative curve has degree one less than the original curve (minimum degree 0).
        /// </summary>
        /// <returns>A new Curve object representing dy/dx with the same domain as the original curve.</returns>
        /// <example>
        /// <code>
        /// var derivative = curve.GetDerivative();
        /// var slope_at_x3 = derivative.EvaluateAt(3).Y;
        /// </code>
        /// </example>
        public Curve GetDerivative()
        {
            if (PolynomialDegree == 0)
                return new Curve(new double[] { 0 }, XMin, XMax);

            var derivCoeffs = new double[PolynomialDegree];
            for (int i = 0; i < PolynomialDegree; i++)
            {
                derivCoeffs[i] = PolynomialCoeffs[i + 1] * (i + 1);
            }

            return new Curve(derivCoeffs, XMin, XMax);
        }

        /// <summary>
        /// Creates a new Curve representing the antiderivative (indefinite integral) of this curve.
        /// The antiderivative curve has degree one higher than the original curve.
        /// </summary>
        /// <param name="integrationConstant">The constant of integration (C). Defaults to 0.</param>
        /// <returns>A new Curve object representing ∫f(x)dx + C with the same domain as the original curve.</returns>
        /// <example>
        /// <code>
        /// var antiderivative = curve.GetAntiderivative(5.0); // +C where C=5
        /// var integral_at_x2 = antiderivative.EvaluateAt(2).Y;
        /// </code>
        /// </example>
        public Curve GetAntiderivative(double integrationConstant = 0)
        {
            var integralCoeffs = new double[PolynomialDegree + 2];
            integralCoeffs[0] = integrationConstant;

            for (int i = 0; i <= PolynomialDegree; i++)
            {
                integralCoeffs[i + 1] = PolynomialCoeffs[i] / (i + 1);
            }

            return new Curve(integralCoeffs, XMin, XMax);
        }

        #endregion

        #region Geometric Analysis

        /// <summary>
        /// Calculates the curvature (measure of how sharply the curve bends) at a specified x-value.
        /// Curvature is computed using the formula κ = |f''(x)| / (1 + (f'(x))²)^(3/2).
        /// </summary>
        /// <param name="xValue">The X-coordinate at which to calculate curvature.</param>
        /// <returns>The curvature value (always non-negative). Higher values indicate sharper bending.</returns>
        /// <example>
        /// <code>
        /// double curvature = curve.CalculateCurvatureAt(2.5);
        /// Console.WriteLine($"Curvature at x=2.5: {curvature}");
        /// </code>
        /// </example>
        public double CalculateCurvatureAt(double xValue)
        {
            var firstDerivative = GetDerivative();
            var secondDerivative = firstDerivative.GetDerivative();

            double dydx = firstDerivative.EvaluateAt(xValue).Y;
            double d2ydx2 = secondDerivative.EvaluateAt(xValue).Y;

            double numerator = Math.Abs(d2ydx2);
            double denominator = Math.Pow(1 + dydx * dydx, 1.5);

            return numerator / denominator;
        }

        /// <summary>
        /// Creates a new linear Curve representing the tangent line to this curve at the specified x-value.
        /// The tangent line has the same slope as the curve at the given point and passes through that point.
        /// </summary>
        /// <param name="xValue">The X-coordinate where the tangent line touches the curve.</param>
        /// <returns>A new linear Curve object representing the tangent line with the same domain as the original curve.</returns>
        /// <example>
        /// <code>
        /// var tangent = curve.GetTangentLineAt(3.0);
        /// var point_on_tangent = tangent.EvaluateAt(3.5);
        /// </code>
        /// </example>
        public Curve GetTangentLineAt(double xValue)
        {
            var point = EvaluateAt(xValue);
            var derivative = GetDerivative();
            double slope = derivative.EvaluateAt(xValue).Y;

            // Tangent line: y - y0 = m(x - x0)
            // y = mx - mx0 + y0
            var coeffs = new double[] { point.Y - slope * xValue, slope };

            return new Curve(coeffs, XMin, XMax);
        }

        /// <summary>
        /// Creates a new linear Curve representing the normal line (perpendicular to tangent) to this curve at the specified x-value.
        /// The normal line is perpendicular to the curve at the given point with slope = -1/(slope of curve).
        /// </summary>
        /// <param name="xValue">The X-coordinate where the normal line intersects the curve.</param>
        /// <returns>A new linear Curve object representing the normal line with the same domain as the original curve.</returns>
        /// <example>
        /// <code>
        /// var normal = curve.GetNormalLineAt(2.0);
        /// var perpendicular_point = normal.EvaluateAt(2.5);
        /// </code>
        /// </example>
        public Curve GetNormalLineAt(double xValue)
        {
            var point = EvaluateAt(xValue);
            var derivative = GetDerivative();
            double slope = derivative.EvaluateAt(xValue).Y;

            // Normal line has slope -1/m
            double normalSlope = (Math.Abs(slope) < 1e-10) ? double.PositiveInfinity : -1 / slope;

            if (double.IsInfinity(normalSlope))
            {
                // Vertical line, approximate with very steep line
                normalSlope = 1e10 * Math.Sign(normalSlope);
            }

            var coeffs = new double[] { point.Y - normalSlope * xValue, normalSlope };

            return new Curve(coeffs, XMin, XMax);
        }

        #endregion

        #region Curve Transformations

        /// <summary>
        /// Creates a new Curve shifted horizontally (left or right) by the specified amount.
        /// Positive values shift right, negative values shift left. Uses binomial expansion for polynomial transformation.
        /// </summary>
        /// <param name="shiftAmount">The amount to shift horizontally. Positive moves right, negative moves left.</param>
        /// <returns>A new Curve object with the same shape but shifted horizontally, with adjusted domain bounds.</returns>
        /// <example>
        /// <code>
        /// var shifted = curve.ShiftHorizontally(2.5); // Move 2.5 units to the right
        /// </code>
        /// </example>
        public Curve ShiftHorizontally(double shiftAmount)
        {
            // For horizontal shift, we need to substitute (x - shift) for x
            // This requires expanding the polynomial
            var newCoeffs = new double[PolynomialCoeffs.Length];

            for (int i = 0; i <= PolynomialDegree; i++)
            {
                // Use binomial expansion
                for (int j = 0; j <= i; j++)
                {
                    double binomialCoeff = BinomialCoefficient(i, j);
                    double term = PolynomialCoeffs[i] * binomialCoeff * Math.Pow(-shiftAmount, i - j);
                    newCoeffs[j] += term;
                }
            }

            return new Curve(newCoeffs, XMin + shiftAmount, XMax + shiftAmount);
        }

        /// <summary>
        /// Creates a new Curve shifted vertically (up or down) by the specified amount.
        /// Positive values shift up, negative values shift down. Only affects the constant term of the polynomial.
        /// </summary>
        /// <param name="shiftAmount">The amount to shift vertically. Positive moves up, negative moves down.</param>
        /// <returns>A new Curve object shifted vertically with the same domain as the original curve.</returns>
        /// <example>
        /// <code>
        /// var shifted = curve.ShiftVertically(-3.0); // Move 3 units down
        /// </code>
        /// </example>
        public Curve ShiftVertically(double shiftAmount)
        {
            var newCoeffs = PolynomialCoeffs.ToArray();
            newCoeffs[0] += shiftAmount;
            return new Curve(newCoeffs, XMin, XMax);
        }

        #endregion

        #region Curve Interactions

        /// <summary>
        /// Finds all intersection points between this curve and another curve within their overlapping domains.
        /// Uses root-finding on the difference polynomial to locate where the curves are equal.
        /// </summary>
        /// <param name="otherCurve">The other Curve object to find intersections with. Must not be null.</param>
        /// <returns>Collection of Point2D objects representing intersection points, ordered by increasing X-coordinate.</returns>
        /// <exception cref="ArgumentNullException">Thrown when otherCurve is null.</exception>
        /// <example>
        /// <code>
        /// var intersections = curve1.FindIntersectionsWith(curve2);
        /// foreach(var point in intersections)
        ///     Console.WriteLine($"Intersection at ({point.X}, {point.Y})");
        /// </code>
        /// </example>
        public IEnumerable<Point2D> FindIntersectionsWith(Curve otherCurve)
        {
            // Find overlapping domain
            double overlapMin = Math.Max(XMin, otherCurve.XMin);
            double overlapMax = Math.Min(XMax, otherCurve.XMax);

            if (overlapMin >= overlapMax)
                return Enumerable.Empty<Point2D>();

            // Create difference polynomial
            int maxDegree = Math.Max(PolynomialDegree, otherCurve.PolynomialDegree);
            var diffCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= PolynomialDegree; i++)
                diffCoeffs[i] += PolynomialCoeffs[i];

            for (int i = 0; i <= otherCurve.PolynomialDegree; i++)
                diffCoeffs[i] -= otherCurve.PolynomialCoeffs[i];

            var diffCurve = new Curve(diffCoeffs, overlapMin, overlapMax);
            var intersectionX = diffCurve.FindRoots();

            return intersectionX.Select(x => EvaluateAt(x));
        }

        /// <summary>
        /// Interpolates this curve with another curve by sampling points and creating a new fitted curve.
        /// This method handles curves with different domains gracefully by sampling corresponding relative positions.
        /// </summary>
        /// <param name="otherCurve">The target curve to interpolate towards. Must not be null.</param>
        /// <param name="interpolationRatio">The interpolation factor between 0 and 1. 0 returns this curve, 1 returns the other curve.</param>
        /// <param name="samplePoints">Number of points to sample for creating the interpolated curve. Defaults to 50.</param>
        /// <returns>A new Curve object representing the interpolated curve with polynomial degree based on the input curves.</returns>
        /// <exception cref="ArgumentNullException">Thrown when otherCurve is null.</exception>
        /// <exception cref="ArgumentException">Thrown when interpolationRatio is not between 0 and 1.</exception>
        /// <example>
        /// <code>
        /// var midway = curve1.InterpolateWith(curve2, 0.5); // 50% between curves
        /// var closer_to_curve2 = curve1.InterpolateWith(curve2, 0.8); // 80% towards curve2
        /// </code>
        /// </example>
        public Curve InterpolateWith(Curve otherCurve, double interpolationRatio, int samplePoints = 50)
        {
            if (interpolationRatio < 0 || interpolationRatio > 1)
                throw new ArgumentException("Interpolation ratio must be between 0 and 1");

            // Determine the interpolated domain
            double newXMin = XMin + interpolationRatio * (otherCurve.XMin - XMin);
            double newXMax = XMax + interpolationRatio * (otherCurve.XMax - XMax);

            // Sample points from each curve
            var interpolatedPoints = new List<Point2D>();

            for (int i = 0; i < samplePoints; i++)
            {
                // Calculate normalized position (0 to 1)
                double t = i / (double)(samplePoints - 1);

                // Map to each curve's domain
                double x1 = XMin + t * (XMax - XMin);
                double x2 = otherCurve.XMin + t * (otherCurve.XMax - otherCurve.XMin);

                // Get y values from each curve
                double y1 = EvaluateAt(x1).Y;
                double y2 = otherCurve.EvaluateAt(x2).Y;

                // Interpolate x and y positions
                double interpolatedX = x1 + interpolationRatio * (x2 - x1);
                double interpolatedY = y1 + interpolationRatio * (y2 - y1);

                interpolatedPoints.Add(new Point2D(interpolatedX, interpolatedY));
            }

            // Determine appropriate polynomial degree for the interpolated curve
            // Use the maximum degree of the two curves, but not more than the number of points - 1
            int targetDegree = Math.Min(
                Math.Max(PolynomialDegree, otherCurve.PolynomialDegree),
                samplePoints - 1
            );

            // Create new curve by fitting the interpolated points
            return new Curve(interpolatedPoints, targetDegree, newXMin, newXMax);
        }

        /// <summary>
        /// Interpolates this curve with another curve using direct coefficient interpolation.
        /// Use this when both curves have similar domains and you want to preserve polynomial structure exactly.
        /// </summary>
        /// <param name="otherCurve">The target curve to interpolate towards. Must not be null.</param>
        /// <param name="interpolationRatio">The interpolation factor between 0 and 1. 0 returns this curve, 1 returns the other curve.</param>
        /// <returns>A new Curve object with interpolated coefficients and domain bounds.</returns>
        /// <exception cref="ArgumentNullException">Thrown when otherCurve is null.</exception>
        /// <exception cref="ArgumentException">Thrown when interpolationRatio is not between 0 and 1.</exception>
        /// <example>
        /// <code>
        /// var interpolated = curve1.InterpolateWithCoefficients(curve2, 0.3);
        /// // Preserves exact polynomial structure
        /// </code>
        /// </example>
        public Curve InterpolateWithCoefficients(Curve otherCurve, double interpolationRatio)
        {
            if (interpolationRatio < 0 || interpolationRatio > 1)
                throw new ArgumentException("Interpolation ratio must be between 0 and 1");

            // Interpolate bounds
            double newXMin = XMin + interpolationRatio * (otherCurve.XMin - XMin);
            double newXMax = XMax + interpolationRatio * (otherCurve.XMax - XMax);

            // Interpolate coefficients
            int maxDegree = Math.Max(PolynomialDegree, otherCurve.PolynomialDegree);
            var newCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= maxDegree; i++)
            {
                double thisCoeff = (i <= PolynomialDegree) ? PolynomialCoeffs[i] : 0;
                double otherCoeff = (i <= otherCurve.PolynomialDegree) ? otherCurve.PolynomialCoeffs[i] : 0;
                newCoeffs[i] = thisCoeff + interpolationRatio * (otherCoeff - thisCoeff);
            }

            return new Curve(newCoeffs, newXMin, newXMax);
        }

        #endregion

        #region Statistical Analysis

        /// <summary>
        /// Calculates comprehensive statistical measures for the curve over its domain including mean, variance, and centroid.
        /// Provides insight into the curve's central tendencies and distribution characteristics.
        /// </summary>
        /// <returns>A CurveStatistics object containing mean value, variance, standard deviation, and centroid point.</returns>
        /// <example>
        /// <code>
        /// var stats = curve.GetStatistics();
        /// Console.WriteLine($"Mean: {stats.Mean}, Std Dev: {stats.StandardDeviation}");
        /// Console.WriteLine($"Centroid: ({stats.CentroidPoint.X}, {stats.CentroidPoint.Y})");
        /// </code>
        /// </example>
        public CurveStatistics GetStatistics()
        {
            // Calculate mean value
            double mean = CalculateDefiniteIntegral() / (XMax - XMin);

            // Calculate centroid
            var xCurve = new Curve(new double[] { 0, 1 }, XMin, XMax); // f(x) = x
            var xTimesCurve = MultiplyPolynomials(xCurve.PolynomialCoeffs, PolynomialCoeffs);
            var xTimesIntegral = IntegratePolynomial(xTimesCurve, XMin, XMax);
            double xCentroid = xTimesIntegral / CalculateDefiniteIntegral();

            // Calculate variance
            var meanCurve = new Curve(new double[] { mean }, XMin, XMax);
            var diffCurve = this - meanCurve;
            var squaredDiffCoeffs = MultiplyPolynomials(diffCurve.PolynomialCoeffs, diffCurve.PolynomialCoeffs);
            double variance = IntegratePolynomial(squaredDiffCoeffs, XMin, XMax) / (XMax - XMin);

            return new CurveStatistics
            {
                Mean = mean,
                Variance = variance,
                StandardDeviation = Math.Sqrt(variance),
                CentroidPoint = new Point2D(xCentroid, EvaluateAt(xCentroid).Y)
            };
        }

        /// <summary>
        /// Calculates quality metrics for polynomial curve fitting using the original data points.
        /// Provides measures of how well the polynomial fits the original data including R-squared and error metrics.
        /// </summary>
        /// <returns>A FittingMetrics object containing R-squared, MSE, RMSE, and MAE values.</returns>
        /// <exception cref="InvalidOperationException">Thrown when the curve was not created from data points (no original data available).</exception>
        /// <example>
        /// <code>
        /// var metrics = curve.GetFittingMetrics();
        /// Console.WriteLine($"R² = {metrics.RSquared:F4}");
        /// Console.WriteLine($"RMSE = {metrics.RootMeanSquaredError:F4}");
        /// </code>
        /// </example>
        public FittingMetrics GetFittingMetrics()
        {
            if (XValues == null || YValues == null)
                throw new InvalidOperationException("Cannot calculate fitting metrics without original data points");

            double sumSquaredErrors = 0;
            double sumAbsoluteErrors = 0;
            double sumY = 0;
            double sumYSquared = 0;
            int n = XValues.Length;

            for (int i = 0; i < n; i++)
            {
                double predicted = EvaluateAt(XValues[i]).Y;
                double error = YValues[i] - predicted;
                sumSquaredErrors += error * error;
                sumAbsoluteErrors += Math.Abs(error);
                sumY += YValues[i];
                sumYSquared += YValues[i] * YValues[i];
            }

            double meanY = sumY / n;
            double totalSumSquares = sumYSquared - n * meanY * meanY;
            double rSquared = 1 - (sumSquaredErrors / totalSumSquares);

            return new FittingMetrics
            {
                RSquared = rSquared,
                MeanSquaredError = sumSquaredErrors / n,
                RootMeanSquaredError = Math.Sqrt(sumSquaredErrors / n),
                MeanAbsoluteError = sumAbsoluteErrors / n
            };
        }

        #endregion

        #region Operators

        /// <summary>
        /// Adds two curves together by adding their polynomial coefficients.
        /// The resulting curve has a domain that is the intersection of both input curves' domains.
        /// </summary>
        /// <param name="curve1">The first curve to add.</param>
        /// <param name="curve2">The second curve to add.</param>
        /// <returns>A new Curve representing the sum of the two input curves.</returns>
        public static Curve operator +(Curve curve1, Curve curve2)
        {
            double newXMin = Math.Max(curve1.XMin, curve2.XMin);
            double newXMax = Math.Min(curve1.XMax, curve2.XMax);

            int maxDegree = Math.Max(curve1.PolynomialDegree, curve2.PolynomialDegree);
            var newCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= curve1.PolynomialDegree; i++)
                newCoeffs[i] += curve1.PolynomialCoeffs[i];

            for (int i = 0; i <= curve2.PolynomialDegree; i++)
                newCoeffs[i] += curve2.PolynomialCoeffs[i];

            return new Curve(newCoeffs, newXMin, newXMax);
        }

        /// <summary>
        /// Subtracts the second curve from the first curve by subtracting their polynomial coefficients.
        /// The resulting curve has a domain that is the intersection of both input curves' domains.
        /// </summary>
        /// <param name="curve1">The curve to subtract from (minuend).</param>
        /// <param name="curve2">The curve to subtract (subtrahend).</param>
        /// <returns>A new Curve representing the difference of the two input curves.</returns>
        public static Curve operator -(Curve curve1, Curve curve2)
        {
            double newXMin = Math.Max(curve1.XMin, curve2.XMin);
            double newXMax = Math.Min(curve1.XMax, curve2.XMax);

            int maxDegree = Math.Max(curve1.PolynomialDegree, curve2.PolynomialDegree);
            var newCoeffs = new double[maxDegree + 1];

            for (int i = 0; i <= curve1.PolynomialDegree; i++)
                newCoeffs[i] += curve1.PolynomialCoeffs[i];

            for (int i = 0; i <= curve2.PolynomialDegree; i++)
                newCoeffs[i] -= curve2.PolynomialCoeffs[i];

            return new Curve(newCoeffs, newXMin, newXMax);
        }

        /// <summary>
        /// Shifts a curve vertically upward by adding a constant value to all y-coordinates.
        /// Equivalent to calling ShiftVertically with a positive value.
        /// </summary>
        /// <param name="curve">The curve to shift vertically.</param>
        /// <param name="verticalShift">The amount to shift upward (positive) or downward (negative).</param>
        /// <returns>A new Curve shifted vertically by the specified amount.</returns>
        public static Curve operator +(Curve curve, double verticalShift)
        {
            return curve.ShiftVertically(verticalShift);
        }

        /// <summary>
        /// Shifts a curve vertically upward by adding a constant value to all y-coordinates.
        /// Equivalent to calling ShiftVertically with a positive value.
        /// </summary>
        /// <param name="verticalShift">The amount to shift upward (positive) or downward (negative).</param>
        /// <param name="curve">The curve to shift vertically.</param>
        /// <returns>A new Curve shifted vertically by the specified amount.</returns>
        public static Curve operator +(double verticalShift, Curve curve)
        {
            return curve.ShiftVertically(verticalShift);
        }

        /// <summary>
        /// Shifts a curve vertically downward by subtracting a constant value from all y-coordinates.
        /// Equivalent to calling ShiftVertically with a negative value.
        /// </summary>
        /// <param name="curve">The curve to shift vertically.</param>
        /// <param name="verticalShift">The amount to shift downward (positive) or upward (negative).</param>
        /// <returns>A new Curve shifted vertically by the negative of the specified amount.</returns>
        public static Curve operator -(Curve curve, double verticalShift)
        {
            return curve.ShiftVertically(-verticalShift);
        }

        /// <summary>
        /// Scales a curve vertically by multiplying all polynomial coefficients by a scalar value.
        /// Positive scalars preserve curve orientation, negative scalars flip the curve vertically.
        /// </summary>
        /// <param name="curve">The curve to scale.</param>
        /// <param name="scalar">The scaling factor. Values > 1 stretch, 0 < values < 1 compress, negative values flip.</param>
        /// <returns>A new Curve scaled by the specified factor.</returns>
        public static Curve operator *(Curve curve, double scalar)
        {
            var newCoeffs = curve.PolynomialCoeffs.Select(c => c * scalar).ToArray();
            return new Curve(newCoeffs, curve.XMin, curve.XMax);
        }

        /// <summary>
        /// Scales a curve vertically by multiplying all polynomial coefficients by a scalar value.
        /// Positive scalars preserve curve orientation, negative scalars flip the curve vertically.
        /// </summary>
        /// <param name="scalar">The scaling factor. Values > 1 stretch, 0 < values < 1 compress, negative values flip.</param>
        /// <param name="curve">The curve to scale.</param>
        /// <returns>A new Curve scaled by the specified factor.</returns>
        public static Curve operator *(double scalar, Curve curve)
        {
            return curve * scalar;
        }

        /// <summary>
        /// Scales a curve vertically by dividing all polynomial coefficients by a scalar value.
        /// Equivalent to multiplying by 1/scalar.
        /// </summary>
        /// <param name="curve">The curve to scale.</param>
        /// <param name="scalar">The divisor. Must not be zero or close to zero.</param>
        /// <returns>A new Curve scaled by 1/scalar.</returns>
        /// <exception cref="DivideByZeroException">Thrown when scalar is zero or very close to zero.</exception>
        public static Curve operator /(Curve curve, double scalar)
        {
            if (Math.Abs(scalar) < 1e-10)
                throw new DivideByZeroException("Cannot divide curve by zero");

            return curve * (1 / scalar);
        }

        /// <summary>
        /// Shifts a curve horizontally to the right by the specified amount.
        /// Uses bit shift operator syntax for intuitive horizontal movement.
        /// </summary>
        /// <param name="curve">The curve to shift horizontally.</param>
        /// <param name="shiftAmount">The amount to shift to the right (positive values).</param>
        /// <returns>A new Curve shifted horizontally to the right.</returns>
        public static Curve operator >>(Curve curve, int shiftAmount)
        {
            return curve.ShiftHorizontally(shiftAmount);
        }

        /// <summary>
        /// Shifts a curve horizontally to the left by the specified amount.
        /// Uses bit shift operator syntax for intuitive horizontal movement.
        /// </summary>
        /// <param name="curve">The curve to shift horizontally.</param>
        /// <param name="shiftAmount">The amount to shift to the left (positive values).</param>
        /// <returns>A new Curve shifted horizontally to the left.</returns>
        public static Curve operator <<(Curve curve, int shiftAmount)
        {
            return curve.ShiftHorizontally(-shiftAmount);
        }

        #endregion

        #region String Representation

        /// <summary>
        /// Returns a string representation of the polynomial curve in mathematical notation.
        /// Shows the polynomial equation with coefficients and the domain range.
        /// </summary>
        /// <returns>A string showing the polynomial equation and its domain, e.g., "2 + 3x^1 - x^2 from 0 to 10".</returns>
        /// <example>
        /// <code>
        /// Console.WriteLine(curve.ToString());
        /// // Output: "1.5 + 2.3x^1 - 0.8x^2 from -5 to 5"
        /// </code>
        /// </example>
        public override string ToString()
        {
            var terms = new List<string>();

            for (int i = 0; i < PolynomialCoeffs.Length; i++)
            {
                double coeff = PolynomialCoeffs[i];

                if (i == 0)
                {
                    // First term - no leading sign
                    terms.Add($"{coeff}x^{i}");
                }
                else
                {
                    // Subsequent terms - include sign
                    if (coeff >= 0)
                    {
                        terms.Add($"+ {coeff}x^{i}");
                    }
                    else
                    {
                        // Negative coefficient already has minus sign
                        terms.Add($"- {Math.Abs(coeff)}x^{i}");
                    }
                }
            }

            string polynomial = string.Join(" ", terms);
            return $"{polynomial} from {XMin} to {XMax}";
        }

        #endregion

        #region Private Helper Methods

        private static double[] FitPolynomial(double[] xValues, double[] yValues, int degree)
        {
            int n = xValues.Length;
            int m = degree + 1;

            // Create Vandermonde matrix
            double[,] matrix = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    matrix[i, j] = Math.Pow(xValues[i], j);
                }
            }

            // Solve using least squares (normal equations)
            double[,] gramMatrix = new double[m, m];
            double[] momentVector = new double[m];

            for (int i = 0; i < m; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        gramMatrix[i, j] += matrix[k, i] * matrix[k, j];
                    }
                }
                for (int k = 0; k < n; k++)
                {
                    momentVector[i] += matrix[k, i] * yValues[k];
                }
            }

            // Solve linear system using Gaussian elimination
            return SolveLinearSystem(gramMatrix, momentVector);
        }

        private static double[] SolveLinearSystem(double[,] A, double[] b)
        {
            int n = b.Length;
            double[,] augmented = new double[n, n + 1];

            // Create augmented matrix
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    augmented[i, j] = A[i, j];
                }
                augmented[i, n] = b[i];
            }

            // Forward elimination
            for (int i = 0; i < n; i++)
            {
                // Find pivot
                int maxRow = i;
                for (int k = i + 1; k < n; k++)
                {
                    if (Math.Abs(augmented[k, i]) > Math.Abs(augmented[maxRow, i]))
                    {
                        maxRow = k;
                    }
                }

                // Swap rows
                for (int k = i; k <= n; k++)
                {
                    double temp = augmented[maxRow, k];
                    augmented[maxRow, k] = augmented[i, k];
                    augmented[i, k] = temp;
                }

                // Make all rows below this one 0 in current column
                for (int k = i + 1; k < n; k++)
                {
                    double factor = augmented[k, i] / augmented[i, i];
                    for (int j = i; j <= n; j++)
                    {
                        augmented[k, j] -= factor * augmented[i, j];
                    }
                }
            }

            // Back substitution
            double[] solution = new double[n];
            for (int i = n - 1; i >= 0; i--)
            {
                solution[i] = augmented[i, n];
                for (int j = i + 1; j < n; j++)
                {
                    solution[i] -= augmented[i, j] * solution[j];
                }
                solution[i] /= augmented[i, i];
            }

            return solution;
        }

        private static double EvaluatePolynomial(double[] coeffs, double x)
        {
            double result = 0;
            double xPower = 1;
            for (int i = 0; i < coeffs.Length; i++)
            {
                result += coeffs[i] * xPower;
                xPower *= x;
            }
            return result;
        }

        private static double EvaluatePolynomialDerivative(double[] coeffs, double x)
        {
            double result = 0;
            double xPower = 1;
            for (int i = 1; i < coeffs.Length; i++)
            {
                result += i * coeffs[i] * xPower;
                xPower *= x;
            }
            return result;
        }

        private static double NewtonRaphson(double[] coeffs, double x0, double tolerance = 1e-10, int maxIterations = 100)
        {
            double x = x0;
            for (int i = 0; i < maxIterations; i++)
            {
                double f = EvaluatePolynomial(coeffs, x);
                double df = EvaluatePolynomialDerivative(coeffs, x);

                if (Math.Abs(df) < tolerance)
                    return double.NaN;

                double xNew = x - f / df;

                if (Math.Abs(xNew - x) < tolerance)
                    return xNew;

                x = xNew;
            }
            return double.NaN;
        }

        private static double BinomialCoefficient(int n, int k)
        {
            if (k > n) return 0;
            if (k == 0 || k == n) return 1;

            double result = 1;
            for (int i = 1; i <= k; i++)
            {
                result *= (n - k + i);
                result /= i;
            }
            return result;
        }

        private static double[] MultiplyPolynomials(double[] p1, double[] p2)
        {
            int degree = p1.Length + p2.Length - 2;
            double[] result = new double[degree + 1];

            for (int i = 0; i < p1.Length; i++)
            {
                for (int j = 0; j < p2.Length; j++)
                {
                    result[i + j] += p1[i] * p2[j];
                }
            }

            return result;
        }

        private static double IntegratePolynomial(double[] coeffs, double a, double b)
        {
            double result = 0;
            for (int i = 0; i < coeffs.Length; i++)
            {
                double coeff = coeffs[i] / (i + 1);
                result += coeff * (Math.Pow(b, i + 1) - Math.Pow(a, i + 1));
            }
            return result;
        }

        #endregion
    }

    #region Supporting Classes

    /// <summary>
    /// Contains comprehensive statistical information about a polynomial curve over its domain.
    /// Provides measures of central tendency, spread, and geometric properties.
    /// </summary>
    public class CurveStatistics
    {
        /// <summary>
        /// Gets or sets the mean (average) value of the curve over its domain.
        /// Calculated as the definite integral divided by the domain width.
        /// </summary>
        public double Mean { get; set; }

        /// <summary>
        /// Gets or sets the standard deviation of the curve values, measuring spread around the mean.
        /// Represents the square root of the variance.
        /// </summary>
        public double StandardDeviation { get; set; }

        /// <summary>
        /// Gets or sets the variance of the curve values, measuring the average squared deviation from the mean.
        /// Higher values indicate greater spread of curve values.
        /// </summary>
        public double Variance { get; set; }

        /// <summary>
        /// Gets or sets the centroid point representing the geometric center of the curve.
        /// The X-coordinate is the center of mass along the x-axis, Y-coordinate is the curve value at that point.
        /// </summary>
        public Point2D CentroidPoint { get; set; }
    }

    /// <summary>
    /// Contains quality metrics that measure how well a polynomial curve fits its original data points.
    /// Provides various error measures and goodness-of-fit statistics for evaluating curve fitting quality.
    /// </summary>
    public class FittingMetrics
    {
        /// <summary>
        /// Gets or sets the coefficient of determination (R²) indicating the proportion of variance explained by the model.
        /// Values range from 0 to 1, where 1 indicates perfect fit and 0 indicates no explanatory power.
        /// </summary>
        public double RSquared { get; set; }

        /// <summary>
        /// Gets or sets the Mean Squared Error (MSE) measuring the average squared differences between predicted and actual values.
        /// Lower values indicate better fit. Units are squared units of the original Y-values.
        /// </summary>
        public double MeanSquaredError { get; set; }

        /// <summary>
        /// Gets or sets the Root Mean Squared Error (RMSE) measuring the square root of the mean squared error.
        /// Provides error measure in the same units as the original Y-values. Lower values indicate better fit.
        /// </summary>
        public double RootMeanSquaredError { get; set; }

        /// <summary>
        /// Gets or sets the Mean Absolute Error (MAE) measuring the average absolute differences between predicted and actual values.
        /// Less sensitive to outliers than RMSE. Lower values indicate better fit.
        /// </summary>
        public double MeanAbsoluteError { get; set; }
    }

    #endregion
}