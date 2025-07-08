# PolyCurve

A lightweight .NET library for polynomial curve generation, manipulation, and analysis. This library can fit polynomials to data points, evaluate them, find roots, compute integrals, locate extrema, and more, all in one place. It uses custom numerical algorithms for polynomial operations and curve fitting.

## Features

- **Polynomial Creation**: Define a polynomial either by direct coefficients or by fitting to (x,y) data.
- **Evaluation**: Evaluate the polynomial at any x-value.
- **Analysis**: Find roots, local/global minima and maxima, inflection points, tangents, normals, etc.
- **Integration & Differentiation**: Automatically compute derivative and antiderivative curves.
- **Domain Operations**: Temporarily change curve bounds with `WithBounds()` and `WithExpandedBounds()`.
- **Advanced Interpolation**: 
  - Standard interpolation and **extrapolation** (ratios beyond [0,1])
  - **Point-constrained interpolation** - create curves that pass through specific points
  - Coefficient-based and sampling-based interpolation methods
- **Operators**: Supports custom operators (+, -, *, /, <<, >>) for polynomial arithmetic and horizontal shifts.
- **Statistical & Fitting Metrics**: Compute mean, variance, standard deviation, centroid, R², RMSE, and more.

## Getting Started 

- Clone/Download this repository or copy the Curve.cs and Point2D.cs into your own .NET project.
- Add a reference to:
  - The PolyCurve namespace/project (if you're keeping it separate).
  - System.Linq for collection operations which you probably already are using.
- Build the solution in Visual Studio (or via dotnet build on the command line, whatever you prefer, you know how these things are).

## Basic usage

```csharp
using PolyCurve;

// Create polynomial from coefficients
double[] coeffs = { -256, 64, 60, -16, 1 }; 
Curve polynomialCurve = new Curve(coeffs, xMin: -3, xMax: 11);

// Basic evaluation
Point2D result = polynomialCurve.EvaluateAt(2.0);
Console.WriteLine($"Value at x=2: {result}");

// Find extrema and roots
Point2D max = polynomialCurve.FindMaximum();
Console.WriteLine($"Global Maximum: {max}");

var roots = polynomialCurve.FindRoots();
Console.WriteLine($"Roots: {string.Join(", ", roots)}");

// Domain operations
var extendedCurve = polynomialCurve.WithBounds(-10, 20);
var moreRoots = extendedCurve.FindRoots();
Console.WriteLine($"Extended roots: {string.Join(", ", moreRoots)}");

// Advanced interpolation
Curve curve1 = new Curve(new double[] { 0, 1 }, 0, 10);      // y = x
Curve curve2 = new Curve(new double[] { 10, -1 }, 0, 10);   // y = 10 - x

// Standard interpolation (50% blend)
Curve blended = curve1.InterpolateWithCoefficients(curve2, 0.5);

// Extrapolation (150% towards curve2)
Curve extrapolated = curve1.InterpolateWithCoefficients(curve2, 1.5);

// Point-constrained interpolation
Point2D targetPoint = new Point2D(3.0, 7.0);
Curve constrained = curve1.InterpolateWithCoefficients(curve2, targetPoint);
```

## Contributing

This library is open source and free for all. Pull requests are welcome! Please make sure any changes are well-documented and (if possible) covered by tests.

## License

Released under the [MIT License](https://github.com/MuratOzc/PolyCurve/blob/master/LICENSE.txt). Feel free to use it for personal or commercial projects.
