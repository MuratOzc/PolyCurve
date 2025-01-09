using CurvesLibrary;

// Example usage

Console.WriteLine("=== Welcome to the Curve Demo ===");

//====================================================
// 1. Creating Curves in Different Ways
//====================================================

// A) Create a curve by directly providing polynomial coefficients:
//    Example polynomial: y = -256 + 64x + 60x^2 - 16x^3 + x^4
double[] polynomialCoeffs = { -256, 64, 60, -16, 1 };
Curve polyCurve = new Curve(polynomialCoeffs, xMin: -3, xMax: 11);
Console.WriteLine("1A) Polynomial curve from coefficients:");
Console.WriteLine($"    {polyCurve}");
Console.WriteLine($"    Domain: [{polyCurve.XMin}, {polyCurve.XMax}]\n");

// B) Create a curve by fitting polynomial to (x,y) data:
//    Let's generate some data from y = sin(x), just as an example
Func<double, double> sinFunc = x => Math.Sin(x);
double[] xData = Enumerable.Range(0, 20).Select(i => i * 0.5).ToArray();  // 0, 0.5, 1.0, 1.5, ...
double[] yData = xData.Select(sinFunc).ToArray();
// Fit a polynomial of degree 5 to these sine data points
Curve fittedCurve = new Curve(xData, yData, polynomialDegree: 5);
Console.WriteLine("1B) Curve from (X,Y) data (sin points), fitted polynomial of degree 5:");
Console.WriteLine($"    Domain: [{fittedCurve.XMin}, {fittedCurve.XMax}]");
Console.WriteLine($"    Coefficients: {fittedCurve}\n");

//====================================================
// 2. Interpolation (Evaluate the curve at an x-value)
//====================================================
double xVal = 2.0;
double yVal = polyCurve.Interpolate(xVal);
Console.WriteLine($"2) Interpolate polyCurve at x={xVal}: y = {yVal}\n");

//====================================================
// 3. InterpolateWith: Merge two curves
//====================================================
// We'll merge polyCurve (1) and fittedCurve (2) with a ratio of 0.5
// (i.e., halfway between them, based on sampling & re-fitting).
Curve mergedCurve = polyCurve.InterpolateWith(fittedCurve, ratio: 0.5);
Console.WriteLine("3) InterpolateWith (50% between polyCurve and fittedCurve)");
Console.WriteLine($"   Merged curve domain: [{mergedCurve.XMin}, {mergedCurve.XMax}]");
Console.WriteLine($"   Merged curve's coeffs: {mergedCurve}\n");

//====================================================
// 4. FindXForY: Solve for x given a y-value
//====================================================
// Suppose we want to find x where polyCurve == 100
var xForY = polyCurve.FindXForY(100);
Console.WriteLine("4) FindXForY(100) for polyCurve");
if (xForY.Any())
{
    foreach (var xFound in xForY)
    {
        Console.WriteLine($"   x={xFound:F4} where curve y=100");
    }
}
else
{
    Console.WriteLine("   No solution found for y=100 in domain.");
}
Console.WriteLine();

//====================================================
// 5. FindMaximum / FindMinimum
//====================================================
var (maxX, maxY) = polyCurve.FindMaximum();
var (minX, minY) = polyCurve.FindMinimum();
Console.WriteLine($"5) polyCurve Maximum at x={maxX:F4}, y={maxY:F4}");
Console.WriteLine($"   polyCurve Minimum at x={minX:F4}, y={minY:F4}\n");

//====================================================
// 6. FindRoots
//====================================================
var roots = polyCurve.FindRoots();
Console.WriteLine("6) polyCurve Roots (x-intercepts):");
if (roots.Any())
{
    foreach (var r in roots)
    {
        Console.WriteLine($"   x={r:F4}");
    }
}
else
{
    Console.WriteLine("   No real roots found in domain.");
}
Console.WriteLine();

//====================================================
// 7. DerivativeCurve and IntegralCurve
//====================================================
var derivative = polyCurve.DerivativeCurve();
Console.WriteLine($"7A) derivative of polyCurve: {derivative}");
var integralCurve = polyCurve.IntegralCurve(c: 10); // add integration constant 10
Console.WriteLine($"7B) integral of polyCurve (+10 constant): {integralCurve}\n");

//====================================================
// 8. FindInflectionPoints
//    Points where the second derivative = 0
//====================================================
var inflections = polyCurve.FindInflectionPoints();
Console.WriteLine("8) Inflection points of polyCurve (x, y):");
if (inflections.Any())
{
    foreach (var (ix, iy) in inflections)
    {
        Console.WriteLine($"   x={ix:F4}, y={iy:F4}");
    }
}
else
{
    Console.WriteLine("   No inflection points found in domain.");
}
Console.WriteLine();

//====================================================
// 9. CalculateCurvature at a point
//====================================================
double curvatureAt2 = polyCurve.CalculateCurvature(2);
Console.WriteLine($"9) Curvature of polyCurve at x=2: {curvatureAt2:F6}\n");

//====================================================
// 10. FindLocalExtrema
//====================================================
var extrema = polyCurve.FindLocalExtrema();
Console.WriteLine("10) Local extrema of polyCurve:");
if (extrema.Any())
{
    foreach (var (lx, ly) in extrema)
    {
        Console.WriteLine($"    x={lx:F4}, y={ly:F4}");
    }
}
else
{
    Console.WriteLine("    No local extrema found in domain.");
}
Console.WriteLine();

//====================================================
// 11. CalculateArcLength (distance along the curve)
//====================================================
double arcLength = polyCurve.CalculateArcLength(-2, 5);
Console.WriteLine($"11) Arc length of polyCurve from x=-2 to x=5: {arcLength:F4}\n");

//====================================================
// 12. CalculateArea vs. CalculateIntegral
//     - CalculateArea: absolute area (|f(x)|)
//     - CalculateIntegral: signed area
//====================================================
double signedArea = polyCurve.CalculateIntegral(-2, 5);
double absoluteArea = polyCurve.CalculateArea(-2, 5);
Console.WriteLine("12) Area under polyCurve from x=-2 to x=5:");
Console.WriteLine($"    Signed area (integral): {signedArea:F4}");
Console.WriteLine($"    Absolute area:          {absoluteArea:F4}\n");

//====================================================
// 13. SplitCurve
//====================================================
double splitPoint = 4.0;
var (leftCurve, rightCurve) = polyCurve.SplitCurve(splitPoint);
Console.WriteLine($"13) Split polyCurve at x={splitPoint}:");
Console.WriteLine($"    Left curve domain:  [{leftCurve.XMin}, {leftCurve.XMax}]");
Console.WriteLine($"    Right curve domain: [{rightCurve.XMin}, {rightCurve.XMax}]\n");

//====================================================
// 14. TangentLine and NormalLine
//====================================================
double xPoint = 2.0;
Curve tangent = polyCurve.TangentLine(xPoint);
Curve normal = polyCurve.NormalLine(xPoint);
Console.WriteLine($"14) Tangent and Normal lines for polyCurve at x={xPoint}:");
Console.WriteLine($"    Tangent: {tangent}");
Console.WriteLine($"    Normal:  {normal}\n");

//====================================================
// 15. ShiftHorizontal
//====================================================
double shiftAmount = 2.0;
var shifted = polyCurve.ShiftHorizontal(shiftAmount);
Console.WriteLine($"15) Shift polyCurve horizontally by +{shiftAmount}:");
Console.WriteLine($"    New domain: [{shifted.XMin}, {shifted.XMax}]");
Console.WriteLine($"    Shifted polynomial: {shifted}\n");

//====================================================
// 16. CalculateStatistics
//====================================================
var stats = polyCurve.CalculateStatistics(50);
Console.WriteLine("16) Statistics (sampled 50 points from polyCurve):");
Console.WriteLine($"    Mean:     {stats.Mean:F4}");
Console.WriteLine($"    Variance: {stats.Variance:F4}");
Console.WriteLine($"    Skewness: {stats.Skewness:F4}");
Console.WriteLine($"    Kurtosis: {stats.Kurtosis:F4}\n");

//====================================================
// 17. CalculateFittingQualityMetrics
//     Only applicable if the curve was constructed with (XValues, YValues).
//====================================================
if (fittedCurve.XValues != null && fittedCurve.YValues != null)
{
    var fitMetrics = fittedCurve.CalculateFittingQualityMetrics();
    Console.WriteLine("17) Fitting Quality Metrics for 'fittedCurve' (fitted to sin data):");
    Console.WriteLine($"    R^2:      {fitMetrics.RSquared:F4}");
    Console.WriteLine($"    MSE:      {fitMetrics.MeanSquaredError:F4}");
    Console.WriteLine($"    RMSE:     {fitMetrics.RootMeanSquaredError:F4}");
    Console.WriteLine($"    MAE:      {fitMetrics.MeanAbsoluteError:F4}\n");
}
else
{
    Console.WriteLine("17) Skipped fitting metrics since no original data for polyCurve.\n");
}

//====================================================
// 18. FindIntersections
//====================================================
// Let's intersect the polynomial curve with a simple line: y = x
double[] lineCoeffs = { 0, 1 }; // line: y = x
Curve lineCurve = new Curve(lineCoeffs, -10, 10);
var intersections = polyCurve.FindIntersections(lineCurve);
Console.WriteLine("18) Intersections between polyCurve and line y=x:");
if (intersections.Any())
{
    foreach (var (ix, iy) in intersections)
    {
        Console.WriteLine($"    Intersection at x={ix:F4}, y={iy:F4}");
    }
}
else
{
    Console.WriteLine("    No intersections found.");
}
Console.WriteLine();

//====================================================
// 19. Operator Overloads (+, -, *, /, etc.)
//====================================================
// We'll demo with smaller polynomials for clarity:
Curve c1 = new Curve(new double[] { 1, 2 }, 0, 5); // y = 1 + 2x
Curve c2 = new Curve(new double[] { -1, 3 }, 0, 5); // y = -1 + 3x

Curve sumCurve = c1 + c2;   // y = (1-1) + (2+3)x = 0 + 5x
Curve diffCurve = c1 - c2;  // y = (1-(-1)) + (2-3)x = 2 + (-1)x
Curve multByScalar = c1 * 2; // y = 2 + 4x
Curve divByScalar = c2 / 2;  // y = -0.5 + 1.5x

Console.WriteLine("19) Operator Overloads:");
Console.WriteLine($"    c1 = {c1}");
Console.WriteLine($"    c2 = {c2}");
Console.WriteLine($"    c1 + c2 = {sumCurve}");
Console.WriteLine($"    c1 - c2 = {diffCurve}");
Console.WriteLine($"    c1 * 2   = {multByScalar}");
Console.WriteLine($"    c2 / 2   = {divByScalar}");

// Scalar addition/subtraction:
Curve addScalar = c1 + 5; // y = (1+5) + 2x = 6 + 2x
Curve subScalar = c2 - 5; // y = (-1-5) + 3x = -6 + 3x
Console.WriteLine($"    c1 + 5   = {addScalar}");
Console.WriteLine($"    c2 - 5   = {subScalar}");

// Shifting operators << and >>:
// c1 >> 2 means shift c1 to the right by +2
// c1 << 2 means shift c1 to the left by 2
Curve shiftedRight = c1 >> 2;
Curve shiftedLeft  = c1 << 2;
Console.WriteLine($"    c1 >> 2  = {shiftedRight}  (Domain: [{shiftedRight.XMin}, {shiftedRight.XMax}])");
Console.WriteLine($"    c1 << 2  = {shiftedLeft}   (Domain: [{shiftedLeft.XMin}, {shiftedLeft.XMax}])\n");

//====================================================
// 20. Done
//====================================================
Console.WriteLine("=== End of Detailed Curve Demo ===");
Console.WriteLine("Press any key to exit.");
Console.ReadKey();