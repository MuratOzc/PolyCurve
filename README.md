
# Curves Library

A lightweight .NET library for polynomial curve generation, manipulation, and analysis. This library can fit polynomials to data points, evaluate them, find roots, compute integrals, locate extrema, and more—all in one place. It uses MathNet.Numerics under the hood in some places and runs custom logic on top of it for easier usage.
## Features

- Polynomial Creation: Define a polynomial either by direct coefficients or by fitting to (x,y) data.
- Evaluation: Interpolate the polynomial at any x within its domain.
- Analysis: Find roots, local/global minima and maxima, inflection points, tangents, normals, etc.
- Integration & Differentiation: Automatically compute derivative and integral curves.
- Operators: Supports custom operators (+, -, *, /, <<, >>) for polynomial arithmetic and horizontal shifts.
- Statistical & Fitting Metrics: Compute mean, variance, skewness, kurtosis, R², MSE, and more.

## Getting Started 

- Clone/Download this repository or copy the Curve.cs (and related files) into your own .NET project.
- Add a reference to:
  - The Curves Library project (if you’re keeping it separate).
  - MathNet.Numerics (via NuGet or local reference) to ensure all numerical operations work properly.
- Build the solution in Visual Studio (or via dotnet build on the command line, whatever you prefer, you know how these things are).

## Basic usage

```csharp
using Curves.Library;

double[] coeffs = { -256, 64, 60, -16, 1 }; 
Curve polynomialCurve = new Curve(coeffs, xMin: -3, xMax: 11);

double result = polynomialCurve.Interpolate(2.0);
Console.WriteLine($"Value at x=2: {result}");

var (maxX, maxY) = polynomialCurve.FindMaximum();
Console.WriteLine($"Global Maximum at x={maxX}, y={maxY}");
```

  
## FAQ (Even though no one asked a question)

#### Does the library handle very large polynomial degrees?

It can, but high-degree polynomials (e.g., >10) may become numerically unstable, and functions like root-finding or extremum detection might be less reliable without additional care.

#### Can I evaluate the curve outside its [XMin, XMax] domain?

By default, calling Interpolate(x) with x outside the domain throws an ArgumentOutOfRangeException. You can extend the domain by creating a new Curve with updated bounds, or handle the exception as needed.

  
## Contributing

This library is open source and free for all. Pull requests are welcome! Please make sure any changes are well-documented and (if possible) covered by tests.
## License

Released under the [MIT License](https://github.com/MuratOzc/CurvesSolution/blob/master/LICENSE.txt). Feel free to use it for personal or commercial projects.
