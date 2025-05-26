using System;
using System.Linq;
using PolyCurve;


Console.WriteLine("=== Comprehensive PolyCurve Library Demonstration ===\n");

//====================================================
// 1. Creating Curves - All Constructor Options
//====================================================
Console.WriteLine("1. CREATING CURVES\n");

// A) From polynomial coefficients
// Polynomial: y = -256 + 64x + 60x^2 - 16x^3 + x^4
double[] coeffs = { -256, 64, 60, -16, 1 };
Curve polynomialCurve = new Curve(coeffs, -3, 11);
Console.WriteLine("1A) Polynomial curve from coefficients:");
Console.WriteLine($"    {polynomialCurve}\n");

// B) From X,Y data points
double[] xData = { 0, 1, 2, 3, 4, 5 };
double[] yData = { 1, 2.5, 5, 8.5, 13, 18.5 };
Curve fittedCurve = new Curve(xData, yData, degree: 2); // Fit quadratic
Console.WriteLine("1B) Curve fitted to data points (degree 2):");
Console.WriteLine($"    {fittedCurve}");

// Show fitting quality
var fitMetrics = fittedCurve.GetFittingMetrics();
Console.WriteLine($"    R² = {fitMetrics.RSquared:F4}, RMSE = {fitMetrics.RootMeanSquaredError:F4}\n");

// C) From Point2D collection
var points = new[] {
            new Point2D(-2, 4),
            new Point2D(-1, 1),
            new Point2D(0, 0),
            new Point2D(1, 1),
            new Point2D(2, 4)
        };
Curve parabolaCurve = new Curve(points, degree: 2);
Console.WriteLine("1C) Curve from Point2D collection:");
Console.WriteLine($"    {parabolaCurve}\n");

//====================================================
// 2. Basic Evaluation
//====================================================
Console.WriteLine("2. BASIC EVALUATION\n");

double x = 5.0;
Point2D point = polynomialCurve.EvaluateAt(x);
Console.WriteLine($"Evaluating polynomial at x = {x}: {point}\n");

//====================================================
// 3. Curve Segmentation
//====================================================
Console.WriteLine("3. CURVE SEGMENTATION\n");

// Split at specific points
var segments1 = polynomialCurve.SplitCurve(0, 4, 8).ToList();
Console.WriteLine("Split at x = 0, 4, 8:");
for (int i = 0; i < segments1.Count; i++)
{
    Console.WriteLine($"    Segment {i + 1}: [{segments1[i].XMin}, {segments1[i].XMax}]");
}

// Split into equal segments
var segments2 = polynomialCurve.SplitCurve(4).ToList();
Console.WriteLine("\nSplit into 4 equal segments:");
for (int i = 0; i < segments2.Count; i++)
{
    Console.WriteLine($"    Segment {i + 1}: [{segments2[i].XMin}, {segments2[i].XMax}]");
}
Console.WriteLine();

//====================================================
// 4. Finding Extrema and Critical Points
//====================================================
Console.WriteLine("4. EXTREMA AND CRITICAL POINTS\n");

Point2D max = polynomialCurve.FindMaximum();
Point2D min = polynomialCurve.FindMinimum();
Console.WriteLine($"Global maximum: {max}");
Console.WriteLine($"Global minimum: {min}");

var localExtrema = polynomialCurve.FindLocalExtrema().ToList();
Console.WriteLine("\nLocal extrema:");
foreach (var extremum in localExtrema)
{
    Console.WriteLine($"    {extremum}");
}

var inflectionPoints = polynomialCurve.FindInflectionPoints().ToList();
Console.WriteLine("\nInflection points:");
foreach (var inflection in inflectionPoints)
{
    Console.WriteLine($"    {inflection}");
}
Console.WriteLine();

//====================================================
// 5. Root Finding
//====================================================
Console.WriteLine("5. ROOT FINDING\n");

var roots = polynomialCurve.FindRoots().ToList();
Console.WriteLine("Roots (x-intercepts):");
foreach (var root in roots)
{
    Console.WriteLine($"    x = {root:F6}");
}

// Find where curve equals specific y-value
double targetY = 100;
var xValues = polynomialCurve.FindXValuesForY(targetY).ToList();
Console.WriteLine($"\nX values where y = {targetY}:");
foreach (var xVal in xValues)
{
    Console.WriteLine($"    x = {xVal:F6}");
}
Console.WriteLine();

//====================================================
// 6. Integration and Area Calculations
//====================================================
Console.WriteLine("6. INTEGRATION AND AREA\n");

// Definite integral over full domain
double integral1 = polynomialCurve.CalculateDefiniteIntegral();
Console.WriteLine($"Integral over full domain [{polynomialCurve.XMin}, {polynomialCurve.XMax}]: {integral1:F2}");

// Definite integral over specific range
double integral2 = polynomialCurve.CalculateDefiniteIntegral(-2, 5);
Console.WriteLine($"Integral from -2 to 5: {integral2:F2}");

// Area under curve (absolute value)
double area = polynomialCurve.CalculateAreaUnderCurve(-2, 5);
Console.WriteLine($"Area under curve from -2 to 5: {area:F2}");

// Arc length
double arcLength = polynomialCurve.CalculateArcLength(-2, 5);
Console.WriteLine($"Arc length from -2 to 5: {arcLength:F2}\n");

//====================================================
// 7. Derivatives and Integrals
//====================================================
Console.WriteLine("7. DERIVATIVES AND INTEGRALS\n");

Curve derivative = polynomialCurve.GetDerivative();
Console.WriteLine($"Derivative: {derivative}");

Curve antiderivative = polynomialCurve.GetAntiderivative(10);
Console.WriteLine($"Antiderivative (C=10): {antiderivative}\n");

//====================================================
// 8. Geometric Analysis
//====================================================
Console.WriteLine("8. GEOMETRIC ANALYSIS\n");

double xPoint = 2.0;
double curvature = polynomialCurve.CalculateCurvatureAt(xPoint);
Console.WriteLine($"Curvature at x = {xPoint}: {curvature:F6}");

Curve tangent = polynomialCurve.GetTangentLineAt(xPoint);
Console.WriteLine($"Tangent line at x = {xPoint}: {tangent}");

Curve normal = polynomialCurve.GetNormalLineAt(xPoint);
Console.WriteLine($"Normal line at x = {xPoint}: {normal}\n");

//====================================================
// 9. Curve Transformations
//====================================================
Console.WriteLine("9. CURVE TRANSFORMATIONS\n");

Curve shiftedH = polynomialCurve.ShiftHorizontally(2);
Console.WriteLine($"Horizontally shifted by +2: {shiftedH}");

Curve shiftedV = polynomialCurve.ShiftVertically(50);
Console.WriteLine($"Vertically shifted by +50: {shiftedV}\n");

//====================================================
// 10. Curve Interactions
//====================================================
Console.WriteLine("10. CURVE INTERACTIONS\n");

// Create another curve for intersection
double[] lineCoeffs = { -100, 20 }; // y = -100 + 20x
Curve lineCurve = new Curve(lineCoeffs, -5, 15);

var intersections = polynomialCurve.FindIntersectionsWith(lineCurve).ToList();
Console.WriteLine($"Intersections with line y = -100 + 20x:");
foreach (var intersection in intersections)
{
    Console.WriteLine($"    {intersection}");
}

// Interpolation between curves
Curve interpolated1 = polynomialCurve.InterpolateWith(lineCurve, 0.3);
Console.WriteLine($"\nInterpolated curve (30% towards line): {interpolated1}");

Curve interpolated2 = polynomialCurve.InterpolateWithCoefficients(lineCurve, 0.3);
Console.WriteLine($"Coefficient interpolation (30% towards line): {interpolated2}\n");

//====================================================
// 11. Statistical Analysis
//====================================================
Console.WriteLine("11. STATISTICAL ANALYSIS\n");

var stats = polynomialCurve.GetStatistics();
Console.WriteLine("Curve statistics:");
Console.WriteLine($"    Mean: {stats.Mean:F4}");
Console.WriteLine($"    Variance: {stats.Variance:F4}");
Console.WriteLine($"    Standard Deviation: {stats.StandardDeviation:F4}");
Console.WriteLine($"    Centroid: {stats.CentroidPoint}\n");

//====================================================
// 12. Operator Overloading Demonstrations
//====================================================
Console.WriteLine("12. OPERATOR OVERLOADING\n");

// Create simple curves for clear demonstration
Curve c1 = new Curve(new double[] { 10, 5, -2 }, 0, 5); // 10 + 5x - 2x²
Curve c2 = new Curve(new double[] { -5, 3, 1 }, 0, 5);  // -5 + 3x + x²

Console.WriteLine($"c1: {c1}");
Console.WriteLine($"c2: {c2}");
Console.WriteLine();

// Arithmetic operations
Curve sum = c1 + c2;
Curve diff = c1 - c2;
Console.WriteLine($"c1 + c2 = {sum}");
Console.WriteLine($"c1 - c2 = {diff}");

// Scalar operations
Curve scaled = c1 * 2;
Curve divided = c2 / 2;
Curve vertUp = c1 + 10;
Curve vertDown = c2 - 5;

Console.WriteLine($"\nc1 * 2 = {scaled}");
Console.WriteLine($"c2 / 2 = {divided}");
Console.WriteLine($"c1 + 10 = {vertUp}");
Console.WriteLine($"c2 - 5 = {vertDown}");

// Shift operators
Curve rightShift = c1 >> 2;
Curve leftShift = c2 << 1;

Console.WriteLine($"\nc1 >> 2 = {rightShift}");
Console.WriteLine($"c2 << 1 = {leftShift}\n");

//====================================================
// 13. Special Cases and Edge Cases
//====================================================
Console.WriteLine("13. SPECIAL CASES\n");

// Constant curve
Curve constant = new Curve(new double[] { 5 }, -10, 10);
Console.WriteLine($"Constant curve: {constant}");
Console.WriteLine($"Derivative of constant: {constant.GetDerivative()}");

// High-degree polynomial from many points
double[] manyX = Enumerable.Range(0, 10).Select(i => (double)i).ToArray();
double[] manyY = manyX.Select(x => Math.Sin(x) * Math.Exp(-x / 5)).ToArray();
Curve complexFit = new Curve(manyX, manyY, degree: -1); // Auto-select degree
Console.WriteLine($"\nComplex curve fitted to sin(x)*exp(-x/5) data:");
Console.WriteLine($"Auto-selected degree: {complexFit.PolynomialDegree}");

//====================================================
// Done
//====================================================
Console.WriteLine("\n=== End of PolyCurve Demonstration ===");
Console.WriteLine("Press any key to exit...");
Console.ReadKey();