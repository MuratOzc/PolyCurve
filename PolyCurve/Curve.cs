using System;
using System.Collections.Generic;
using System.Linq;

namespace PolyCurve
{
    /// <summary>
    /// Represents a polynomial curve segment with various mathematical operations and analysis methods.
    /// </summary>
    public class Curve
    {
        public double[]? XValues { get; }
        public double[]? YValues { get; }
        public double XMin { get; }
        public double XMax { get; }
        public int PolynomialDegree { get; }
        public double[] PolynomialCoeffs { get; }

        #region Constructors

        /// <summary>
        /// Creates a polynomial curve by fitting data points with specified degree within given bounds.
        /// </summary>
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
        /// Creates a polynomial curve by fitting a collection of points with specified degree within given bounds.
        /// </summary>
        public Curve(IEnumerable<Point2D> dataPoints, int degree = -1, double? xMin = null, double? xMax = null)
            : this(dataPoints.Select(p => p.X).ToArray(),
                   dataPoints.Select(p => p.Y).ToArray(),
                   degree, xMin, xMax)
        {
        }

        /// <summary>
        /// Creates a polynomial curve directly from coefficients with specified bounds.
        /// </summary>
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
        /// Evaluates the polynomial at the specified x-value and returns the corresponding point.
        /// </summary>
        public Point2D EvaluateAt(double xValue)
        {
            double y = EvaluatePolynomial(PolynomialCoeffs, xValue);
            return new Point2D(xValue, y);
        }

        #endregion

        #region Curve Segmentation

        /// <summary>
        /// Splits the curve at the specified x-values, creating multiple curve segments.
        /// </summary>
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
        /// Splits the curve into equally-sized segments.
        /// </summary>
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
        /// Finds the maximum point on the curve within its domain.
        /// </summary>
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
        /// Finds the minimum point on the curve within its domain.
        /// </summary>
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
        /// Finds all local extrema (minima and maxima) within the curve's domain.
        /// </summary>
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
        /// Finds the inflection points of the curve within its domain.
        /// </summary>
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
        /// Finds the roots (x-intercepts) of the curve within its domain.
        /// </summary>
        public IEnumerable<double> FindRoots()
        {
            return FindXValuesForY(0);
        }

        /// <summary>
        /// Finds the x-values corresponding to a given y-value within the curve's domain.
        /// </summary>
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
        /// Calculates the definite integral of the curve over its domain.
        /// </summary>
        public double CalculateDefiniteIntegral()
        {
            return CalculateDefiniteIntegral(XMin, XMax);
        }

        /// <summary>
        /// Calculates the definite integral of the curve between specified bounds.
        /// </summary>
        public double CalculateDefiniteIntegral(double xStart, double xEnd)
        {
            var antiderivative = GetAntiderivative();
            double endValue = antiderivative.EvaluateAt(xEnd).Y;
            double startValue = antiderivative.EvaluateAt(xStart).Y;
            return endValue - startValue;
        }

        /// <summary>
        /// Calculates the area under the curve over its domain.
        /// </summary>
        public double CalculateAreaUnderCurve()
        {
            return CalculateAreaUnderCurve(XMin, XMax);
        }

        /// <summary>
        /// Calculates the area under the curve between specified bounds.
        /// </summary>
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
        /// Calculates the arc length of the curve over its domain.
        /// </summary>
        public double CalculateArcLength()
        {
            return CalculateArcLength(XMin, XMax);
        }

        /// <summary>
        /// Calculates the arc length of the curve between specified bounds.
        /// </summary>
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
        /// Creates a new Curve representing the derivative of this curve.
        /// </summary>
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
        /// Creates a new Curve representing the antiderivative of this curve.
        /// </summary>
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
        /// Calculates the curvature of the curve at a given x-value.
        /// </summary>
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
        /// Creates a new Curve representing the tangent line at the specified x-value.
        /// </summary>
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
        /// Creates a new Curve representing the normal line at the specified x-value.
        /// </summary>
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
        /// Creates a new Curve shifted horizontally by the specified amount.
        /// </summary>
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
        /// Creates a new Curve shifted vertically by the specified amount.
        /// </summary>
        public Curve ShiftVertically(double shiftAmount)
        {
            var newCoeffs = PolynomialCoeffs.ToArray();
            newCoeffs[0] += shiftAmount;
            return new Curve(newCoeffs, XMin, XMax);
        }

        #endregion

        #region Curve Interactions

        /// <summary>
        /// Finds the intersection points between this curve and another curve within overlapping domains.
        /// </summary>
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
        /// This method handles curves with different domains more gracefully than coefficient interpolation.
        /// </summary>
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
        /// Interpolates this curve with another curve using coefficient interpolation.
        /// Use this when both curves have similar domains and you want to preserve polynomial structure.
        /// </summary>
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
        /// Calculates various statistical measures for the curve over its domain.
        /// </summary>
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
        /// Calculates quality metrics for the polynomial fitting using original data points.
        /// </summary>
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

        // Curve operations
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

        // Vertical shifting
        public static Curve operator +(Curve curve, double verticalShift)
        {
            return curve.ShiftVertically(verticalShift);
        }

        public static Curve operator +(double verticalShift, Curve curve)
        {
            return curve.ShiftVertically(verticalShift);
        }

        public static Curve operator -(Curve curve, double verticalShift)
        {
            return curve.ShiftVertically(-verticalShift);
        }

        // Scaling
        public static Curve operator *(Curve curve, double scalar)
        {
            var newCoeffs = curve.PolynomialCoeffs.Select(c => c * scalar).ToArray();
            return new Curve(newCoeffs, curve.XMin, curve.XMax);
        }

        public static Curve operator *(double scalar, Curve curve)
        {
            return curve * scalar;
        }

        public static Curve operator /(Curve curve, double scalar)
        {
            if (Math.Abs(scalar) < 1e-10)
                throw new DivideByZeroException("Cannot divide curve by zero");

            return curve * (1 / scalar);
        }

        // Horizontal shifting
        public static Curve operator >>(Curve curve, int shiftAmount)
        {
            return curve.ShiftHorizontally(shiftAmount);
        }

        public static Curve operator <<(Curve curve, int shiftAmount)
        {
            return curve.ShiftHorizontally(-shiftAmount);
        }

        #endregion

        #region String Representation

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
    /// Contains statistical information about a polynomial curve.
    /// </summary>
    public class CurveStatistics
    {
        public double Mean { get; set; }
        public double StandardDeviation { get; set; }
        public double Variance { get; set; }
        public Point2D CentroidPoint { get; set; }
    }

    /// <summary>
    /// Contains metrics about the quality of polynomial curve fitting.
    /// </summary>
    public class FittingMetrics
    {
        public double RSquared { get; set; }
        public double MeanSquaredError { get; set; }
        public double RootMeanSquaredError { get; set; }
        public double MeanAbsoluteError { get; set; }
    }

    #endregion
}