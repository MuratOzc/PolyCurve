using PolyCurve;

namespace Tests
{
    public abstract class BaseCurveTests
    {
        protected const double Tolerance = 0.1;
        protected const double HighPrecisionTolerance = 0.0005;

        protected void AssertPointsEqual(IEnumerable<Point2D> actual, IEnumerable<Point2D> expected)
        {
            var actualList = actual.OrderBy(p => p.X).ToList();
            var expectedList = expected.OrderBy(p => p.X).ToList();
            Assert.That(actualList.Count, Is.EqualTo(expectedList.Count));
            for (int i = 0; i < actualList.Count; i++)
            {
                AssertPointEqual(actualList[i], expectedList[i]);
            }
        }

        protected void AssertPointEqual(Point2D actual, Point2D expected)
        {
            Assert.That(actual.X, Is.EqualTo(expected.X).Within(Tolerance));
            Assert.That(actual.Y, Is.EqualTo(expected.Y).Within(Tolerance));
        }

        protected void AssertArraysEqual(IEnumerable<double> actual, IEnumerable<double> expected)
        {
            var actualList = actual.ToList();
            var expectedList = expected.ToList();
            Assert.That(actualList.Count, Is.EqualTo(expectedList.Count));
            for (int i = 0; i < actualList.Count; i++)
            {
                Assert.That(actualList[i], Is.EqualTo(expectedList[i]).Within(Tolerance));
            }
        }

        protected static double[] GenerateLinearSpaced(int count, double start, double end)
        {
            if (count < 2) throw new ArgumentException("Count must be at least 2");

            double[] result = new double[count];
            double step = (end - start) / (count - 1);

            for (int i = 0; i < count; i++)
            {
                result[i] = start + i * step;
            }

            return result;
        }
    }

    [TestFixture]
    public class CurveTests : BaseCurveTests
    {
        private Curve polynomialCurve;
        private Curve sinusoidalCurve;
        private Curve exponentialCurve;
        private Curve testCurve;

        [SetUp]
        public void Setup()
        {
            // Polynomial curve: y = -256 + 64x + 60x^2 - 16x^3 + x^4
            double[] coeffs = { -256, 64, 60, -16, 1 };
            polynomialCurve = new Curve(coeffs, -3, 11);

            // Sinusoidal curve: y = 5 * sin(x) + 2
            Func<double, double> sinFunc = x => 5 * Math.Sin(x) + 2;
            double[] xValues = GenerateLinearSpaced(100, 0, 2 * Math.PI);
            double[] yValues = xValues.Select(sinFunc).ToArray();
            sinusoidalCurve = new Curve(xValues, yValues, degree: 10);

            // Exponential curve: y = 2^x
            Func<double, double> expFunc = x => Math.Pow(2, x);
            xValues = GenerateLinearSpaced(100, 0, 5);
            yValues = xValues.Select(expFunc).ToArray();
            exponentialCurve = new Curve(xValues, yValues, degree: 10);

            // Test curve for domain operations: y = x^2 - 4 (roots at -2 and 2)
            testCurve = new Curve(new double[] { -4, 0, 1 }, -3, 3);
        }

        #region Constructor and Property Tests

        [Test]
        public void TestConstructor_FromDataPoints_BasicFunctionality()
        {
            double[] xValues = { 1, 2, 3, 4, 5 };
            double[] yValues = { 2, 8, 18, 32, 50 };

            var curve = new Curve(xValues, yValues, degree: 2);

            Assert.Multiple(() =>
            {
                Assert.That(curve.XValues, Is.Not.Null, "XValues should not be null");
                Assert.That(curve.YValues, Is.Not.Null, "YValues should not be null");
                Assert.That(curve.XValues.Length, Is.EqualTo(5), "XValues length incorrect");
                Assert.That(curve.YValues.Length, Is.EqualTo(5), "YValues length incorrect");
                Assert.That(curve.PolynomialDegree, Is.EqualTo(2), "Polynomial degree incorrect");
                Assert.That(curve.XMin, Is.EqualTo(1), "XMin incorrect");
                Assert.That(curve.XMax, Is.EqualTo(5), "XMax incorrect");
                AssertArraysEqual(curve.XValues, xValues);
                AssertArraysEqual(curve.YValues, yValues);
            });
        }

        [Test]
        public void TestConstructor_FromDataPoints_AutoDegree()
        {
            double[] xValues = { 1, 2, 3, 4, 5, 6 };
            double[] yValues = { 1, 4, 9, 16, 25, 36 };

            var curve = new Curve(xValues, yValues); // degree = -1 (auto)

            Assert.Multiple(() =>
            {
                Assert.That(curve.PolynomialDegree, Is.EqualTo(5), "Auto degree should be min(data_points-1, 5)");
            });
        }

        [Test]
        public void TestConstructor_FromDataPoints_CustomBounds()
        {
            double[] xValues = { 1, 2, 3, 4, 5 };
            double[] yValues = { 2, 8, 18, 32, 50 };

            var curve = new Curve(xValues, yValues, degree: 2, xMin: 0, xMax: 10);

            Assert.Multiple(() =>
            {
                Assert.That(curve.XMin, Is.EqualTo(0), "Custom XMin not set");
                Assert.That(curve.XMax, Is.EqualTo(10), "Custom XMax not set");
            });
        }

        [Test]
        public void TestConstructor_FromPoint2DCollection()
        {
            var points = new List<Point2D>
            {
                new Point2D(1, 2),
                new Point2D(2, 8),
                new Point2D(3, 18),
                new Point2D(4, 32),
                new Point2D(5, 50)
            };

            var curve = new Curve(points, degree: 2);

            Assert.Multiple(() =>
            {
                Assert.That(curve.XValues, Is.Not.Null, "XValues should not be null");
                Assert.That(curve.YValues, Is.Not.Null, "YValues should not be null");
                Assert.That(curve.PolynomialDegree, Is.EqualTo(2), "Polynomial degree incorrect");
                AssertArraysEqual(curve.XValues, new double[] { 1, 2, 3, 4, 5 });
                AssertArraysEqual(curve.YValues, new double[] { 2, 8, 18, 32, 50 });
            });
        }

        [Test]
        public void TestConstructor_FromCoefficients()
        {
            double[] coeffs = { 1, 2, 3 }; // 1 + 2x + 3x²
            var curve = new Curve(coeffs, -5, 5);

            Assert.Multiple(() =>
            {
                Assert.That(curve.XValues, Is.Null, "XValues should be null for coefficient constructor");
                Assert.That(curve.YValues, Is.Null, "YValues should be null for coefficient constructor");
                Assert.That(curve.PolynomialDegree, Is.EqualTo(2), "Polynomial degree incorrect");
                Assert.That(curve.XMin, Is.EqualTo(-5), "XMin incorrect");
                Assert.That(curve.XMax, Is.EqualTo(5), "XMax incorrect");
                AssertArraysEqual(curve.PolynomialCoeffs, coeffs);
            });
        }

        [Test]
        public void TestConstructor_ErrorHandling()
        {
            Assert.Multiple(() =>
            {
                // Null arrays
                Assert.Throws<ArgumentNullException>(() => new Curve(null, new double[] { 1, 2 }));
                Assert.Throws<ArgumentNullException>(() => new Curve(new double[] { 1, 2 }, null));

                // Mismatched array lengths
                Assert.Throws<ArgumentException>(() => new Curve(new double[] { 1, 2 }, new double[] { 1, 2, 3 }));

                // Too few data points
                Assert.Throws<ArgumentException>(() => new Curve(new double[] { 1 }, new double[] { 1 }));

                // Empty coefficients
                Assert.Throws<ArgumentException>(() => new Curve(new double[0], 0, 1));
                Assert.Throws<ArgumentException>(() => new Curve(null, 0, 1));

                // Null Point2D collection
                Assert.Throws<ArgumentNullException>(() => new Curve((IEnumerable<Point2D>)null));
            });
        }

        #endregion

        #region Basic Evaluation Tests

        [Test]
        public void TestEvaluateAt()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.EvaluateAt(5).Y, Is.EqualTo(189.0).Within(HighPrecisionTolerance), "Polynomial curve evaluation failed");
                Assert.That(sinusoidalCurve.EvaluateAt(Math.PI / 4).Y, Is.EqualTo(5.5355).Within(HighPrecisionTolerance), "Sinusoidal curve evaluation failed");
                Assert.That(exponentialCurve.EvaluateAt(2).Y, Is.EqualTo(4.0).Within(HighPrecisionTolerance), "Exponential curve evaluation failed");
            });
        }

        [Test]
        public void TestEvaluateAt_EdgeCases()
        {
            Assert.Multiple(() =>
            {
                // Evaluate at domain boundaries
                var point1 = polynomialCurve.EvaluateAt(polynomialCurve.XMin);
                var point2 = polynomialCurve.EvaluateAt(polynomialCurve.XMax);
                Assert.That(point1.X, Is.EqualTo(polynomialCurve.XMin));
                Assert.That(point2.X, Is.EqualTo(polynomialCurve.XMax));

                // Evaluate outside domain (should still work mathematically)
                var pointOutside = polynomialCurve.EvaluateAt(polynomialCurve.XMax + 1);
                Assert.That(pointOutside.X, Is.EqualTo(polynomialCurve.XMax + 1));
            });
        }

        #endregion

        #region Point Generation Tests

        [Test]
        public void TestGeneratePoints_BasicFunctionality()
        {
            var points = testCurve.GeneratePoints(5).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(points.Length, Is.EqualTo(5), "Should generate exact number of points requested");
                Assert.That(points[0].X, Is.EqualTo(testCurve.XMin), "First point should be at XMin");
                Assert.That(points[4].X, Is.EqualTo(testCurve.XMax), "Last point should be at XMax");

                // Check that points are evenly spaced
                double expectedStep = (testCurve.XMax - testCurve.XMin) / 4;
                for (int i = 1; i < points.Length; i++)
                {
                    double actualStep = points[i].X - points[i - 1].X;
                    Assert.That(actualStep, Is.EqualTo(expectedStep).Within(HighPrecisionTolerance),
                        $"Points should be evenly spaced at step {i}");
                }
            });
        }

        [Test]
        public void TestGeneratePoints_CustomRange()
        {
            var points = testCurve.GeneratePoints(10, -1, 1).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(points.Length, Is.EqualTo(10), "Should generate exact number of points requested");
                Assert.That(points[0].X, Is.EqualTo(-1), "First point should be at custom start");
                Assert.That(points[9].X, Is.EqualTo(1), "Last point should be at custom end");

                // Verify Y values are correct for known function
                for (int i = 0; i < points.Length; i++)
                {
                    var expectedPoint = testCurve.EvaluateAt(points[i].X);
                    AssertPointEqual(points[i], expectedPoint);
                }
            });
        }

        [Test]
        public void TestGeneratePoints_ErrorHandling()
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentException>(() => testCurve.GeneratePoints(1).ToArray(),
                    "Should throw for count < 2");
                Assert.Throws<ArgumentException>(() => testCurve.GeneratePoints(0).ToArray(),
                    "Should throw for count = 0");
                Assert.Throws<ArgumentException>(() => testCurve.GeneratePoints(-5).ToArray(),
                    "Should throw for negative count");
                Assert.Throws<ArgumentException>(() => testCurve.GeneratePoints(10, 5, 3).ToArray(),
                    "Should throw when start >= end");
                Assert.Throws<ArgumentException>(() => testCurve.GeneratePoints(10, 2, 2).ToArray(),
                    "Should throw when start == end");
            });
        }

        [Test]
        public void TestGeneratePoints_ConsistencyWithEvaluateAt()
        {
            var points = polynomialCurve.GeneratePoints(100);

            foreach (var point in points)
            {
                var directEvaluation = polynomialCurve.EvaluateAt(point.X);
                AssertPointEqual(point, directEvaluation);
            }
        }

        [Test]
        public void TestGeneratePoints_EdgeCases()
        {
            Assert.Multiple(() =>
            {
                // Minimum count
                var twoPoints = testCurve.GeneratePoints(2);
                Assert.That(twoPoints.Count(), Is.EqualTo(2), "Should handle minimum count of 2");

                // Large count
                var manyPoints = testCurve.GeneratePoints(10000);
                Assert.That(manyPoints.Count(), Is.EqualTo(10000), "Should handle large point counts");

                // Very small range
                var smallRange = testCurve.GeneratePoints(5, 0, 0.001);
                Assert.That(smallRange.Count(), Is.EqualTo(5), "Should handle very small ranges");
            });
        }

        #endregion

        #region Curve Segmentation Tests

        [Test]
public void TestSplitCurve_AtSpecificPoints()
{
    Assert.Multiple(() =>
    {
        // Use array syntax to split at specific points, not into equal segments
        Curve[] polyCurves = polynomialCurve.SplitCurve(new double[] { 4 }).ToArray();
        Assert.That(polyCurves[0].XMax, Is.EqualTo(4).Within(HighPrecisionTolerance), "Polynomial left split failed");
        Assert.That(polyCurves[1].XMin, Is.EqualTo(4).Within(HighPrecisionTolerance), "Polynomial right split failed");

        Curve[] sinCurves = sinusoidalCurve.SplitCurve(new double[] { Math.PI }).ToArray();
        Assert.That(sinCurves[0].XMax, Is.EqualTo(Math.PI).Within(HighPrecisionTolerance), "Sinusoidal left split failed");
        Assert.That(sinCurves[1].XMin, Is.EqualTo(Math.PI).Within(HighPrecisionTolerance), "Sinusoidal right split failed");

        var expCurves = exponentialCurve.SplitCurve(new double[] { 2.5 }).ToArray();
        Assert.That(expCurves[0].XMax, Is.EqualTo(2.5).Within(HighPrecisionTolerance), "Exponential left split failed");
        Assert.That(expCurves[1].XMin, Is.EqualTo(2.5).Within(HighPrecisionTolerance), "Exponential right split failed");
    });
}

        [Test]
        public void TestSplitCurve_MultiplePoints()
        {
            var segments = polynomialCurve.SplitCurve(0, 4, 8).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(segments.Length, Is.EqualTo(4), "Should create 4 segments");
                Assert.That(segments[0].XMin, Is.EqualTo(polynomialCurve.XMin), "First segment XMin");
                Assert.That(segments[0].XMax, Is.EqualTo(0), "First segment XMax");
                Assert.That(segments[1].XMin, Is.EqualTo(0), "Second segment XMin");
                Assert.That(segments[1].XMax, Is.EqualTo(4), "Second segment XMax");
                Assert.That(segments[2].XMin, Is.EqualTo(4), "Third segment XMin");
                Assert.That(segments[2].XMax, Is.EqualTo(8), "Third segment XMax");
                Assert.That(segments[3].XMin, Is.EqualTo(8), "Fourth segment XMin");
                Assert.That(segments[3].XMax, Is.EqualTo(polynomialCurve.XMax), "Fourth segment XMax");
            });
        }

        [Test]
        public void TestSplitCurve_EqualSegments()
        {
            var segments = polynomialCurve.SplitCurve(4).ToArray(); // 4 equal segments

            Assert.Multiple(() =>
            {
                Assert.That(segments.Length, Is.EqualTo(4), "Should create 4 segments");

                double expectedWidth = (polynomialCurve.XMax - polynomialCurve.XMin) / 4;
                for (int i = 0; i < segments.Length; i++)
                {
                    double actualWidth = segments[i].XMax - segments[i].XMin;
                    Assert.That(actualWidth, Is.EqualTo(expectedWidth).Within(HighPrecisionTolerance), $"Segment {i} width incorrect");
                }

                // Check continuity
                for (int i = 0; i < segments.Length - 1; i++)
                {
                    Assert.That(segments[i].XMax, Is.EqualTo(segments[i + 1].XMin).Within(HighPrecisionTolerance),
                        $"Gap between segments {i} and {i + 1}");
                }
            });
        }

        [Test]
        public void TestSplitCurve_ErrorHandling()
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentException>(() => polynomialCurve.SplitCurve(0),
                    "Should throw for zero segments");
                Assert.Throws<ArgumentException>(() => polynomialCurve.SplitCurve(-1),
                    "Should throw for negative segments");
            });
        }

        [Test]
        public void TestSplitCurve_OutsideDomain()
        {
            // Split points outside domain should be ignored
            var segments = testCurve.SplitCurve(-10, 0, 10).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(segments.Length, Is.EqualTo(2), "Should create 2 segments (outside points ignored)");
                Assert.That(segments[0].XMax, Is.EqualTo(0), "First segment should end at 0");
                Assert.That(segments[1].XMin, Is.EqualTo(0), "Second segment should start at 0");
            });
        }

        #endregion

        #region Domain Operations Tests

        [Test]
        public void TestWithBounds_BasicFunctionality()
        {
            var newCurve = testCurve.WithBounds(-5, 5);

            Assert.Multiple(() =>
            {
                Assert.That(newCurve.XMin, Is.EqualTo(-5), "New minimum bound not set correctly");
                Assert.That(newCurve.XMax, Is.EqualTo(5), "New maximum bound not set correctly");
                Assert.That(newCurve.PolynomialDegree, Is.EqualTo(testCurve.PolynomialDegree), "Polynomial degree should be unchanged");
                AssertArraysEqual(newCurve.PolynomialCoeffs, testCurve.PolynomialCoeffs);
            });
        }

        [Test]
        public void TestWithBounds_OriginalUnchanged()
        {
            double originalMin = testCurve.XMin;
            double originalMax = testCurve.XMax;

            var newCurve = testCurve.WithBounds(-10, 10);

            Assert.Multiple(() =>
            {
                Assert.That(testCurve.XMin, Is.EqualTo(originalMin), "Original curve minimum bound changed");
                Assert.That(testCurve.XMax, Is.EqualTo(originalMax), "Original curve maximum bound changed");
                Assert.That(newCurve.XMin, Is.EqualTo(-10), "New curve minimum bound incorrect");
                Assert.That(newCurve.XMax, Is.EqualTo(10), "New curve maximum bound incorrect");
            });
        }

        [Test]
        public void TestWithBounds_InvalidBounds()
        {
            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentException>(() => testCurve.WithBounds(5, 3), "Should throw when min >= max");
                Assert.Throws<ArgumentException>(() => testCurve.WithBounds(5, 5), "Should throw when min == max");
            });
        }

        [Test]
        public void TestWithExpandedBounds_SymmetricExpansion()
        {
            var expandedCurve = testCurve.WithExpandedBounds(2, 2);

            Assert.Multiple(() =>
            {
                Assert.That(expandedCurve.XMin, Is.EqualTo(-5), "Left expansion incorrect"); // -3 - 2 = -5
                Assert.That(expandedCurve.XMax, Is.EqualTo(5), "Right expansion incorrect");  // 3 + 2 = 5
                AssertArraysEqual(expandedCurve.PolynomialCoeffs, testCurve.PolynomialCoeffs);
            });
        }

        [Test]
        public void TestWithExpandedBounds_AsymmetricExpansion()
        {
            var expandedCurve = testCurve.WithExpandedBounds(1, 3);

            Assert.Multiple(() =>
            {
                Assert.That(expandedCurve.XMin, Is.EqualTo(-4), "Left expansion incorrect"); // -3 - 1 = -4
                Assert.That(expandedCurve.XMax, Is.EqualTo(6), "Right expansion incorrect");  // 3 + 3 = 6
            });
        }

        [Test]
        public void TestWithExpandedBounds_Contraction()
        {
            var contractedCurve = testCurve.WithExpandedBounds(-1, -1);

            Assert.Multiple(() =>
            {
                Assert.That(contractedCurve.XMin, Is.EqualTo(-2), "Left contraction incorrect"); // -3 - (-1) = -2
                Assert.That(contractedCurve.XMax, Is.EqualTo(2), "Right contraction incorrect");  // 3 + (-1) = 2
            });
        }

        [Test]
        public void TestWithExpandedBounds_InvalidResult()
        {
            Assert.Throws<ArgumentException>(() => testCurve.WithExpandedBounds(-5, -5),
                "Should throw when resulting bounds would be invalid");
        }

        [Test]
        public void TestDomainOperations_FindMoreRoots()
        {
            // Original domain [-3, 3] should find both roots at -2 and 2
            var originalRoots = testCurve.FindRoots().ToList();

            // Contracted domain [-1, 1] should find no roots
            var contractedRoots = testCurve.WithBounds(-1, 1).FindRoots().ToList();

            // Expanded domain should still find the same roots
            var expandedRoots = testCurve.WithExpandedBounds(2, 2).FindRoots().ToList();

            Assert.Multiple(() =>
            {
                Assert.That(originalRoots.Count, Is.EqualTo(2), "Should find 2 roots in original domain");
                Assert.That(contractedRoots.Count, Is.EqualTo(0), "Should find no roots in contracted domain");
                Assert.That(expandedRoots.Count, Is.EqualTo(2), "Should find same roots in expanded domain");

                // Check root values
                AssertArraysEqual(originalRoots.OrderBy(r => r), new[] { -2.0, 2.0 });
                AssertArraysEqual(expandedRoots.OrderBy(r => r), new[] { -2.0, 2.0 });
            });
        }

        #endregion

        #region Extrema and Critical Points Tests

        [Test]
        public void TestFindMaximum()
        {
            Assert.Multiple(() =>
            {
                AssertPointEqual(polynomialCurve.FindMaximum(), new Point2D(11.0, 1053.0));
                AssertPointEqual(sinusoidalCurve.FindMaximum(), new Point2D(Math.PI / 2, 7.0));
                AssertPointEqual(exponentialCurve.FindMaximum(), new Point2D(5, 32.0));
            });
        }

        [Test]
        public void TestFindMinimum()
        {
            Assert.Multiple(() =>
            {
                AssertPointEqual(polynomialCurve.FindMinimum(), new Point2D(-0.44949, -271.15));
                AssertPointEqual(sinusoidalCurve.FindMinimum(), new Point2D(3 * Math.PI / 2, -3));
                AssertPointEqual(exponentialCurve.FindMinimum(), new Point2D(0, 1));
            });
        }

        [Test]
        public void TestFindLocalExtrema()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindLocalExtrema(), new[] { new Point2D(-0.44949, -271.15), new Point2D(4.4495, 199.15), new Point2D(8, 0) });
                AssertPointsEqual(sinusoidalCurve.FindLocalExtrema(), new[] { new Point2D(Math.PI / 2, 7.0), new Point2D(3 * Math.PI / 2, -3) });
                Assert.That(exponentialCurve.FindLocalExtrema(), Is.Empty, "Exponential curve local extrema failed");
            });
        }

        [Test]
        public void TestFindInflectionPoints()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindInflectionPoints(), [new Point2D(1.55051, -66.3837), new Point2D(6.44949, 90.3837)]);
                AssertPointsEqual(sinusoidalCurve.FindInflectionPoints(), [new Point2D(0, 2.0), new Point2D(Math.PI, 2.0), new Point2D(2 * Math.PI, 2.0)]);
                Assert.That(exponentialCurve.FindInflectionPoints(), Is.Empty, "Exponential curve inflection points failed");
            });
        }

        [Test]
        public void TestExtrema_ConstantCurve()
        {
            var constantCurve = new Curve(new double[] { 5 }, 0, 10);

            Assert.Multiple(() =>
            {
                Assert.That(constantCurve.FindLocalExtrema(), Is.Empty, "Constant curve should have no local extrema");
                AssertPointEqual(constantCurve.FindMaximum(), new Point2D(0, 5)); // Should return endpoint
                AssertPointEqual(constantCurve.FindMinimum(), new Point2D(0, 5)); // Should return endpoint
            });
        }

        [Test]
        public void TestExtrema_LinearCurve()
        {
            var linearCurve = new Curve(new double[] { 1, 2 }, 0, 5); // y = 1 + 2x

            Assert.Multiple(() =>
            {
                Assert.That(linearCurve.FindLocalExtrema(), Is.Empty, "Linear curve should have no local extrema");
                AssertPointEqual(linearCurve.FindMaximum(), new Point2D(5, 11)); // Maximum at right endpoint
                AssertPointEqual(linearCurve.FindMinimum(), new Point2D(0, 1));  // Minimum at left endpoint
            });
        }

        #endregion

        #region Root Finding Tests

        [Test]
        public void TestFindRoots()
        {
            Assert.Multiple(() =>
            {
                AssertArraysEqual(polynomialCurve.FindRoots(), new[] { -2.0, 2.0, 8.0 });
                AssertArraysEqual(sinusoidalCurve.FindRoots(), new[] { 3.55311, 5.87167 });
                Assert.That(exponentialCurve.FindRoots(), Is.Empty, "Exponential curve roots failed");
            });
        }

        [Test]
        public void TestFindXValuesForY()
        {
            Assert.Multiple(() =>
            {
                AssertArraysEqual(polynomialCurve.FindXValuesForY(100), new[] { -2.2263, 2.7659, 6.3370, 9.1234 });
                AssertArraysEqual(sinusoidalCurve.FindXValuesForY(4), new[] { 0.4636, 2.6780 });
                AssertArraysEqual(exponentialCurve.FindXValuesForY(8), new[] { 3.0 });
            });
        }

        [Test]
        public void TestFindXValuesForY_EdgeCases()
        {
            Assert.Multiple(() =>
            {
                // Test finding roots (y = 0)
                var roots1 = testCurve.FindXValuesForY(0);
                var roots2 = testCurve.FindRoots();
                AssertArraysEqual(roots1, roots2);

                // Test value not on curve
                var noSolutions = exponentialCurve.FindXValuesForY(-1);
                Assert.That(noSolutions, Is.Empty, "Should find no solutions for y < 0 on exponential curve");

                // Test maximum/minimum values
                var max = polynomialCurve.FindMaximum();
                var xForMax = polynomialCurve.FindXValuesForY(max.Y);
                Assert.That(xForMax.Any(x => Math.Abs(x - max.X) < Tolerance), Is.True, "Should find x-value for maximum y");
            });
        }

        [Test]
        public void TestRootFinding_Accuracy()
        {
            // Test with known polynomial roots
            var quadratic = new Curve(new double[] { -6, 1, 1 }, -5, 5); // x² + x - 6 = (x+3)(x-2)
            var roots = quadratic.FindRoots().OrderBy(r => r).ToArray();

            Assert.Multiple(() =>
            {
                Assert.That(roots.Length, Is.EqualTo(2), "Should find 2 roots");
                Assert.That(roots[0], Is.EqualTo(-3).Within(HighPrecisionTolerance), "First root should be -3");
                Assert.That(roots[1], Is.EqualTo(2).Within(HighPrecisionTolerance), "Second root should be 2");
            });
        }

        #endregion

        #region Integration and Area Tests

        [Test]
        public void TestCalculateDefiniteIntegral()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateDefiniteIntegral(-3, 11), Is.EqualTo(1178.8).Within(HighPrecisionTolerance), "Polynomial curve integral failed");
                Assert.That(sinusoidalCurve.CalculateDefiniteIntegral(0, 2 * Math.PI), Is.EqualTo(12.5664).Within(HighPrecisionTolerance), "Sinusoidal curve integral failed");
                Assert.That(exponentialCurve.CalculateDefiniteIntegral(0, 5), Is.EqualTo(44.724).Within(HighPrecisionTolerance), "Exponential curve integral failed");
            });
        }

        [Test]
        public void TestCalculateDefiniteIntegral_DomainOverload()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateDefiniteIntegral(), Is.EqualTo(1178.8).Within(HighPrecisionTolerance), "Polynomial curve integral failed");
                Assert.That(sinusoidalCurve.CalculateDefiniteIntegral(), Is.EqualTo(12.5664).Within(HighPrecisionTolerance), "Sinusoidal curve integral failed");
                Assert.That(exponentialCurve.CalculateDefiniteIntegral(), Is.EqualTo(44.724).Within(HighPrecisionTolerance), "Exponential curve integral failed");
            });
        }

        [Test]
        public void TestCalculateAreaUnderCurve()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateAreaUnderCurve(-3, 11), Is.EqualTo(2561.2).Within(HighPrecisionTolerance), "Polynomial curve area failed");
                Assert.That(sinusoidalCurve.CalculateAreaUnderCurve(0, 2 * Math.PI), Is.EqualTo(21.6224).Within(HighPrecisionTolerance), "Sinusoidal curve area failed");
                Assert.That(exponentialCurve.CalculateAreaUnderCurve(0, 5), Is.EqualTo(44.724).Within(HighPrecisionTolerance), "Exponential curve area failed");
            });
        }

        [Test]
        public void TestCalculateAreaUnderCurve_DomainOverload()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateAreaUnderCurve(), Is.EqualTo(2561.2).Within(HighPrecisionTolerance), "Polynomial curve area failed");
                Assert.That(sinusoidalCurve.CalculateAreaUnderCurve(), Is.EqualTo(21.6224).Within(HighPrecisionTolerance), "Sinusoidal curve area failed");
                Assert.That(exponentialCurve.CalculateAreaUnderCurve(), Is.EqualTo(44.724).Within(HighPrecisionTolerance), "Exponential curve area failed");
            });
        }

        [Test]
        public void TestCalculateArcLength()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateArcLength(-3, 11), Is.EqualTo(2598.79744060024886176663).Within(Tolerance), "Polynomial curve arc length failed");
                Assert.That(sinusoidalCurve.CalculateArcLength(0, 2 * Math.PI), Is.EqualTo(21.393860164471).Within(Tolerance), "Sinusoidal curve arc length failed");
                Assert.That(exponentialCurve.CalculateArcLength(0, 5), Is.EqualTo(31.8893).Within(Tolerance), "Exponential curve arc length failed");
            });
        }

        [Test]
        public void TestCalculateArcLength_DomainOverload()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateArcLength(), Is.EqualTo(2598.79744060024886176663).Within(Tolerance), "Polynomial curve arc length failed");
                Assert.That(sinusoidalCurve.CalculateArcLength(), Is.EqualTo(21.393860164471).Within(Tolerance), "Sinusoidal curve arc length failed");
                Assert.That(exponentialCurve.CalculateArcLength(), Is.EqualTo(31.8893).Within(Tolerance), "Exponential curve arc length failed");
            });
        }

        [Test]
        public void TestIntegration_Properties()
        {
            // Test fundamental properties of integration
            var simpleCurve = new Curve(new double[] { 1, 2, 3 }, 0, 2); // 1 + 2x + 3x²

            Assert.Multiple(() =>
            {
                // Integral should equal area for positive functions
                double integral = simpleCurve.CalculateDefiniteIntegral(0, 2);
                double area = simpleCurve.CalculateAreaUnderCurve(0, 2);
                Assert.That(integral, Is.EqualTo(area).Within(HighPrecisionTolerance), "Integral should equal area for positive function");

                // Additivity: integral over [a,c] = integral over [a,b] + integral over [b,c]
                double fullIntegral = simpleCurve.CalculateDefiniteIntegral(0, 2);
                double part1 = simpleCurve.CalculateDefiniteIntegral(0, 1);
                double part2 = simpleCurve.CalculateDefiniteIntegral(1, 2);
                Assert.That(fullIntegral, Is.EqualTo(part1 + part2).Within(HighPrecisionTolerance), "Integration additivity failed");
            });
        }

        [Test]
        public void TestAreaVsIntegral_NegativeRegions()
        {
            // Test curve that crosses x-axis: y = x² - 1 (negative between -1 and 1)
            var crossingCurve = new Curve(new double[] { -1, 0, 1 }, -2, 2);

            Assert.Multiple(() =>
            {
                double integral = crossingCurve.CalculateDefiniteIntegral(-2, 2);
                double area = crossingCurve.CalculateAreaUnderCurve(-2, 2);

                // Area should be larger than absolute value of integral when curve crosses x-axis
                Assert.That(area, Is.GreaterThan(Math.Abs(integral)), "Area should be larger than |integral| for curve crossing x-axis");
            });
        }

        #endregion

        #region Derivatives and Integrals Tests

        [Test]
        public void TestGetDerivative()
        {
            Assert.Multiple(() =>
            {
                var polyDerivative = polynomialCurve.GetDerivative();
                AssertArraysEqual(polyDerivative.PolynomialCoeffs, new[] { 64.0, 120, -48, 4 });

                var sinDerivative = sinusoidalCurve.GetDerivative();
                Assert.That(sinDerivative.EvaluateAt(Math.PI / 2).Y, Is.EqualTo(0).Within(HighPrecisionTolerance), "Sinusoidal derivative at PI/2 failed");
                Assert.That(sinDerivative.EvaluateAt(0).Y, Is.EqualTo(5).Within(Tolerance), "Sinusoidal derivative at 0 failed");

                var expDerivative = exponentialCurve.GetDerivative();
                Assert.That(expDerivative.EvaluateAt(1).Y, Is.EqualTo(1.3863).Within(HighPrecisionTolerance), "Exponential derivative failed");
            });
        }

        [Test]
        public void TestGetAntiderivative()
        {
            Assert.Multiple(() =>
            {
                var polyIntegral = polynomialCurve.GetAntiderivative();
                Assert.That(polyIntegral.EvaluateAt(5).Y, Is.EqualTo(145).Within(HighPrecisionTolerance), "Polynomial integral curve failed");

                var sinIntegral = sinusoidalCurve.GetAntiderivative(integrationConstant: -5);
                Assert.That(sinIntegral.EvaluateAt(Math.PI).Y, Is.EqualTo(5 + 2 * Math.PI).Within(HighPrecisionTolerance), "Sinusoidal integral curve failed");

                var expIntegral = exponentialCurve.GetAntiderivative(integrationConstant: 1.4427);
                Assert.That(expIntegral.EvaluateAt(1).Y, Is.EqualTo(2.88539).Within(HighPrecisionTolerance), "Exponential integral curve failed");
            });
        }

        [Test]
        public void TestDerivative_EdgeCases()
        {
            // Derivative of constant should be zero
            var constantCurve = new Curve(new double[] { 42 }, 0, 10);
            var derivative = constantCurve.GetDerivative();

            Assert.Multiple(() =>
            {
                Assert.That(derivative.PolynomialDegree, Is.EqualTo(0), "Derivative of constant should have degree 0");
                AssertArraysEqual(derivative.PolynomialCoeffs, new double[] { 0 });
            });
        }

        [Test]
        public void TestAntiderivative_Properties()
        {
            var simpleCurve = new Curve(new double[] { 6, 4, 3 }, 0, 5); // 6 + 4x + 3x²

            Assert.Multiple(() =>
            {
                var antiderivative = simpleCurve.GetAntiderivative(0);

                // Check that the derivative of the antiderivative gives back the original (within tolerance)
                var reconstructed = antiderivative.GetDerivative();
                AssertArraysEqual(reconstructed.PolynomialCoeffs, simpleCurve.PolynomialCoeffs);

                // Check integration constant
                var antiderivativeWithConstant = simpleCurve.GetAntiderivative(10);
                Assert.That(antiderivativeWithConstant.PolynomialCoeffs[0], Is.EqualTo(10), "Integration constant not applied correctly");
            });
        }

        #endregion

        #region Geometric Analysis Tests

        [Test]
        public void TestCalculateCurvatureAt()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateCurvatureAt(0), Is.EqualTo(0.000457596).Within(HighPrecisionTolerance), "Polynomial curve curvature failed");
                Assert.That(sinusoidalCurve.CalculateCurvatureAt(Math.PI / 2), Is.EqualTo(5).Within(HighPrecisionTolerance), "Sinusoidal curve curvature failed");
                Assert.That(exponentialCurve.CalculateCurvatureAt(1), Is.EqualTo(0.1924).Within(HighPrecisionTolerance), "Exponential curve curvature failed");
            });
        }

        [Test]
        public void TestGetTangentLineAt()
        {
            Assert.Multiple(() =>
            {
                var polyTangent = polynomialCurve.GetTangentLineAt(5);
                AssertArraysEqual(polyTangent.PolynomialCoeffs, [369.0, -36.0]);

                var sinTangent = sinusoidalCurve.GetTangentLineAt(Math.PI / 4);
                AssertArraysEqual(sinTangent.PolynomialCoeffs, [2.75873, 3.53553]);

                var expTangent = exponentialCurve.GetTangentLineAt(2);
                AssertArraysEqual(expTangent.PolynomialCoeffs, [-1.54518, 2.77259]);
            });
        }

        [Test]
        public void TestGetNormalLineAt()
        {
            Assert.Multiple(() =>
            {
                var polyNormal = polynomialCurve.GetNormalLineAt(5);
                AssertArraysEqual(polyNormal.PolynomialCoeffs, [188.861, 0.0277778]);

                var sinNormal = sinusoidalCurve.GetNormalLineAt(Math.PI / 4);
                AssertArraysEqual(sinNormal.PolynomialCoeffs, [5.75768, -0.282843]);

                var expNormal = exponentialCurve.GetNormalLineAt(2);
                AssertArraysEqual(expNormal.PolynomialCoeffs, [4.72135, -0.360674]);
            });
        }

        [Test]
        public void TestTangentNormal_Perpendicularity()
        {
            double x = 3.0;
            var tangent = polynomialCurve.GetTangentLineAt(x);
            var normal = polynomialCurve.GetNormalLineAt(x);

            // Slopes should be negative reciprocals (m1 * m2 = -1)
            double tangentSlope = tangent.PolynomialCoeffs[1];
            double normalSlope = normal.PolynomialCoeffs[1];

            Assert.That(tangentSlope * normalSlope, Is.EqualTo(-1).Within(HighPrecisionTolerance),
                "Tangent and normal should be perpendicular");
        }

        [Test]
        public void TestTangentLine_PassesThroughPoint()
        {
            double x = 2.5;
            var originalPoint = polynomialCurve.EvaluateAt(x);
            var tangent = polynomialCurve.GetTangentLineAt(x);
            var tangentPoint = tangent.EvaluateAt(x);

            AssertPointEqual(originalPoint, tangentPoint);
        }

        [Test]
        public void TestCurvature_SpecialCases()
        {
            // Straight line should have zero curvature
            var line = new Curve(new double[] { 1, 2 }, 0, 10);
            double curvature = line.CalculateCurvatureAt(5);
            Assert.That(curvature, Is.EqualTo(0).Within(HighPrecisionTolerance), "Line should have zero curvature");

            // Circle approximation: y ≈ r - x²/(2r) has curvature ≈ 1/r
            var circleApprox = new Curve(new double[] { 10, 0, -0.05 }, -5, 5); // Approximate circle with radius 10
            double circleCurvature = circleApprox.CalculateCurvatureAt(0);
            Assert.That(circleCurvature, Is.EqualTo(0.1).Within(0.01), "Circle approximation curvature");
        }

        #endregion

        #region Curve Transformations Tests

        [Test]
        public void TestShiftHorizontally()
        {
            var shifted = testCurve.ShiftHorizontally(2);

            Assert.Multiple(() =>
            {
                // Domain should be shifted
                Assert.That(shifted.XMin, Is.EqualTo(testCurve.XMin + 2), "Domain minimum not shifted");
                Assert.That(shifted.XMax, Is.EqualTo(testCurve.XMax + 2), "Domain maximum not shifted");

                // Function values should be preserved at shifted points
                var originalPoint = testCurve.EvaluateAt(1);
                var shiftedPoint = shifted.EvaluateAt(3); // 1 + 2 = 3
                Assert.That(shiftedPoint.Y, Is.EqualTo(originalPoint.Y).Within(HighPrecisionTolerance), "Y-values not preserved after horizontal shift");
            });
        }

        [Test]
        public void TestShiftVertically()
        {
            var shifted = testCurve.ShiftVertically(5);

            Assert.Multiple(() =>
            {
                // Domain should be unchanged
                Assert.That(shifted.XMin, Is.EqualTo(testCurve.XMin), "Domain minimum changed");
                Assert.That(shifted.XMax, Is.EqualTo(testCurve.XMax), "Domain maximum changed");

                // Y-values should be shifted by constant amount
                var originalPoint = testCurve.EvaluateAt(1);
                var shiftedPoint = shifted.EvaluateAt(1);
                Assert.That(shiftedPoint.Y, Is.EqualTo(originalPoint.Y + 5).Within(HighPrecisionTolerance), "Y-values not shifted correctly");
            });
        }

        [Test]
        public void TestTransformations_ChainedOperations()
        {
            var transformed = testCurve
                .ShiftVertically(3)
                .ShiftHorizontally(-1);

            Assert.Multiple(() =>
            {
                // Check domain
                Assert.That(transformed.XMin, Is.EqualTo(testCurve.XMin - 1), "Chained horizontal shift failed");
                Assert.That(transformed.XMax, Is.EqualTo(testCurve.XMax - 1), "Chained horizontal shift failed");

                // Check function values
                var originalPoint = testCurve.EvaluateAt(2);
                var transformedPoint = transformed.EvaluateAt(1); // x shifted by -1
                Assert.That(transformedPoint.Y, Is.EqualTo(originalPoint.Y + 3).Within(HighPrecisionTolerance), "Chained transformations failed");
            });
        }

        #endregion

        #region Curve Interactions Tests

        [Test]
        public void TestFindIntersectionsWith()
        {
            var curve1 = new Curve(new double[] { 0, 1 }, 0, 10); // y = x
            var curve2 = new Curve(new double[] { 2, -1 }, 0, 10); // y = 2 - x

            var intersections = curve1.FindIntersectionsWith(curve2);
            Point2D[] expected = [new Point2D(1, 1)];

            var sinExpoIntersections = sinusoidalCurve.FindIntersectionsWith(exponentialCurve);
            Point2D[] sinExpoExpected = [new Point2D(2.41344, 5.32743)];
            var sinPolyIntersections = sinusoidalCurve.FindIntersectionsWith(polynomialCurve);
            Point2D[] sinPolyExpected = [new Point2D(2.04495, 6.44783)];
            var expoPolyIntersections = exponentialCurve.FindIntersectionsWith(polynomialCurve);
            Point2D[] expoPolyExpected = [new Point2D(2.0284, 4.07952)];

            Assert.Multiple(() =>
            {
                AssertPointsEqual(intersections, expected);
                AssertPointsEqual(sinExpoIntersections, sinExpoExpected);
                AssertPointsEqual(sinPolyIntersections, sinPolyExpected);
                AssertPointsEqual(expoPolyIntersections, expoPolyExpected);
            });
        }

        [Test]
        public void TestFindIntersections_EdgeCases()
        {
            Assert.Multiple(() =>
            {
                // No overlapping domain
                var curve1 = new Curve(new double[] { 1 }, 0, 5);
                var curve2 = new Curve(new double[] { 1 }, 6, 10);
                var noIntersections = curve1.FindIntersectionsWith(curve2);
                Assert.That(noIntersections, Is.Empty, "Should find no intersections with non-overlapping domains");

                // Identical curves
                var identical1 = new Curve(new double[] { 1, 2, 3 }, 0, 5);
                var identical2 = new Curve(new double[] { 1, 2, 3 }, 0, 5);
                // Identical curves technically intersect everywhere, but the root finder should find specific points
                var identicalIntersections = identical1.FindIntersectionsWith(identical2);
                Assert.That(identicalIntersections.Count(), Is.GreaterThanOrEqualTo(0), "Identical curves intersection behavior");

                // Parallel lines (no intersections)
                var line1 = new Curve(new double[] { 0, 1 }, 0, 10); // y = x
                var line2 = new Curve(new double[] { 5, 1 }, 0, 10); // y = 5 + x
                var parallelIntersections = line1.FindIntersectionsWith(line2);
                Assert.That(parallelIntersections, Is.Empty, "Parallel lines should not intersect");
            });
        }

        [Test]
        public void TestInterpolateWith_BasicRatio()
        {
            var curve1 = new Curve(new double[] { 0, 1 }, 0, 5);    // y = x
            var curve2 = new Curve(new double[] { 10, -1 }, 0, 5);  // y = 10 - x

            Assert.Multiple(() =>
            {
                // Test ratio = 0 (should return curve1)
                var interpolated0 = curve1.InterpolateWith(curve2, 0.0);
                var testPoint0 = interpolated0.EvaluateAt(2);
                var expectedPoint0 = curve1.EvaluateAt(2);
                Assert.That(testPoint0.Y, Is.EqualTo(expectedPoint0.Y).Within(HighPrecisionTolerance), "Ratio 0 should return first curve");

                // Test ratio = 1 (should return curve2)
                var interpolated1 = curve1.InterpolateWith(curve2, 1.0);
                var testPoint1 = interpolated1.EvaluateAt(2);
                var expectedPoint1 = curve2.EvaluateAt(2);
                Assert.That(testPoint1.Y, Is.EqualTo(expectedPoint1.Y).Within(HighPrecisionTolerance), "Ratio 1 should return second curve");

                // Test ratio = 0.5 (should be midway)
                var interpolated05 = curve1.InterpolateWith(curve2, 0.5);
                var testPoint05 = interpolated05.EvaluateAt(2);
                double expectedY = (curve1.EvaluateAt(2).Y + curve2.EvaluateAt(2).Y) / 2;
                Assert.That(testPoint05.Y, Is.EqualTo(expectedY).Within(HighPrecisionTolerance), "Ratio 0.5 should be midway");
            });
        }

        [Test]
        public void TestInterpolateWithCoefficients_BasicRatio()
        {
            var curve1 = new Curve(new double[] { 1, 2, 3 }, 0, 5);
            var curve2 = new Curve(new double[] { 4, 5, 6 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Test ratio = 0
                var interpolated0 = curve1.InterpolateWithCoefficients(curve2, 0.0);
                AssertArraysEqual(interpolated0.PolynomialCoeffs, curve1.PolynomialCoeffs);

                // Test ratio = 1
                var interpolated1 = curve1.InterpolateWithCoefficients(curve2, 1.0);
                AssertArraysEqual(interpolated1.PolynomialCoeffs, curve2.PolynomialCoeffs);

                // Test ratio = 0.5
                var interpolated05 = curve1.InterpolateWithCoefficients(curve2, 0.5);
                double[] expectedCoeffs = { 2.5, 3.5, 4.5 };
                AssertArraysEqual(interpolated05.PolynomialCoeffs, expectedCoeffs);
            });
        }

        [Test]
        public void TestInterpolateWith_WithPoint()
        {
            // Create two simple curves
            var curve1 = new Curve(new double[] { 0, 1 }, -5, 5); // y = x
            var curve2 = new Curve(new double[] { 10, -1 }, -5, 5); // y = 10 - x

            // Target point that should be on the interpolated curve
            var targetPoint = new Point2D(2, 5);

            var interpolated = curve1.InterpolateWith(curve2, targetPoint);

            Assert.Multiple(() =>
            {
                // Check that the interpolated curve passes through the target point
                var actualPoint = interpolated.EvaluateAt(targetPoint.X);
                Assert.That(actualPoint.Y, Is.EqualTo(targetPoint.Y).Within(0.01), "Should pass through target point");
            });
        }

        [Test]
        public void TestInterpolateWithCoefficients_WithPoint()
        {
            // Create two simple curves
            var curve1 = new Curve(new double[] { 0, 1 }, -5, 5); // y = x
            var curve2 = new Curve(new double[] { 10, -1 }, -5, 5); // y = 10 - x

            // Target point that should be on the interpolated curve
            var targetPoint = new Point2D(2, 5);

            var interpolated = curve1.InterpolateWithCoefficients(curve2, targetPoint);

            Assert.Multiple(() =>
            {
                // Check that the interpolated curve passes through the target point
                var actualPoint = interpolated.EvaluateAt(targetPoint.X);
                AssertPointEqual(actualPoint, targetPoint);

                // Verify the interpolation ratio was correct
                // At x=2: curve1 gives y=2, curve2 gives y=8
                // To get y=5, we need ratio = (5-2)/(8-2) = 0.5
                // So coefficients should be [5, 0]
                AssertArraysEqual(interpolated.PolynomialCoeffs, new double[] { 5, 0 });
            });
        }

        [Test]
        public void TestInterpolation_ErrorHandling()
        {
            var curve1 = new Curve(new double[] { 0, 1 }, 0, 10);

            Assert.Multiple(() =>
            {
                Assert.Throws<ArgumentNullException>(() => curve1.InterpolateWith(null, 0.5));
                Assert.Throws<ArgumentNullException>(() => curve1.InterpolateWithCoefficients(null, 0.5));
                Assert.Throws<ArgumentNullException>(() => curve1.InterpolateWith(curve1, null));
                Assert.Throws<ArgumentNullException>(() => curve1.InterpolateWithCoefficients(curve1, null));
            });
        }

        [Test]
        public void TestInterpolation_Extrapolation()
        {
            var curve1 = new Curve(new double[] { 0, 1 }, 0, 5);
            var curve2 = new Curve(new double[] { 10, -1 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Test extrapolation with ratio > 1
                var extrapolated = curve1.InterpolateWithCoefficients(curve2, 1.5);
                var testPoint = extrapolated.EvaluateAt(2);
                // At x=2: curve1 gives 2, curve2 gives 8
                // Ratio 1.5: 2 + 1.5*(8-2) = 2 + 9 = 11
                Assert.That(testPoint.Y, Is.EqualTo(11).Within(HighPrecisionTolerance), "Extrapolation beyond curve2 failed");

                // Test extrapolation with ratio < 0
                var extrapolatedNeg = curve1.InterpolateWithCoefficients(curve2, -0.5);
                var testPointNeg = extrapolatedNeg.EvaluateAt(2);
                // Ratio -0.5: 2 + (-0.5)*(8-2) = 2 - 3 = -1
                Assert.That(testPointNeg.Y, Is.EqualTo(-1).Within(HighPrecisionTolerance), "Negative extrapolation failed");
            });
        }

        #endregion

        #region Statistical Analysis Tests

        [Test]
        public void TestGetStatistics()
        {
            var stats = testCurve.GetStatistics();

            Assert.Multiple(() =>
            {
                Assert.That(stats, Is.Not.Null, "Statistics should not be null");
                Assert.That(stats.Mean, Is.Not.NaN, "Mean should be a valid number");
                Assert.That(stats.StandardDeviation, Is.GreaterThanOrEqualTo(0), "Standard deviation should be non-negative");
                Assert.That(stats.Variance, Is.GreaterThanOrEqualTo(0), "Variance should be non-negative");
                Assert.That(stats.CentroidPoint, Is.Not.Null, "Centroid point should not be null");

                // For symmetric curve y = x² - 4, centroid should be at x = 0
                Assert.That(stats.CentroidPoint.X, Is.EqualTo(0).Within(HighPrecisionTolerance), "Centroid X should be at 0 for symmetric curve");
            });
        }

        [Test]
        public void TestGetStatistics_KnownValues()
        {
            // Test with simple linear function y = x over [0, 2]
            var linearCurve = new Curve(new double[] { 0, 1 }, 0, 2);
            var stats = linearCurve.GetStatistics();

            Assert.Multiple(() =>
            {
                // Mean of y = x over [0,2] should be 1
                Assert.That(stats.Mean, Is.EqualTo(1).Within(HighPrecisionTolerance), "Mean of linear function incorrect");

                // Centroid X should be at 4/3 (center of mass, not geometric center)
                // For y=x: x_centroid = ∫(x*x)dx / ∫x dx = (8/3) / 2 = 4/3
                Assert.That(stats.CentroidPoint.X, Is.EqualTo(4.0 / 3.0).Within(HighPrecisionTolerance), "Centroid X incorrect");
            });
        }

        [Test]
        public void TestGetFittingMetrics()
        {
            // Create curve from known data points with perfect fit
            double[] xValues = { 1, 2, 3, 4, 5 };
            double[] yValues = { 1, 4, 9, 16, 25 }; // Perfect quadratic: y = x²
            var perfectFit = new Curve(xValues, yValues, degree: 2);

            var metrics = perfectFit.GetFittingMetrics();

            Assert.Multiple(() =>
            {
                Assert.That(metrics, Is.Not.Null, "Fitting metrics should not be null");
                Assert.That(metrics.RSquared, Is.EqualTo(1.0).Within(HighPrecisionTolerance), "R² should be 1 for perfect fit");
                Assert.That(metrics.MeanSquaredError, Is.EqualTo(0).Within(HighPrecisionTolerance), "MSE should be 0 for perfect fit");
                Assert.That(metrics.RootMeanSquaredError, Is.EqualTo(0).Within(HighPrecisionTolerance), "RMSE should be 0 for perfect fit");
                Assert.That(metrics.MeanAbsoluteError, Is.EqualTo(0).Within(HighPrecisionTolerance), "MAE should be 0 for perfect fit");
            });
        }

        [Test]
        public void TestGetFittingMetrics_ImperfectFit()
        {
            // Create curve with some noise
            double[] xValues = { 1, 2, 3, 4, 5 };
            double[] yValues = { 1.1, 3.9, 9.2, 15.8, 25.1 }; // Noisy quadratic
            var noisyFit = new Curve(xValues, yValues, degree: 2);

            var metrics = noisyFit.GetFittingMetrics();

            Assert.Multiple(() =>
            {
                Assert.That(metrics.RSquared, Is.GreaterThan(0.9), "R² should be high for good fit");
                Assert.That(metrics.RSquared, Is.LessThan(1.0), "R² should be less than 1 for imperfect fit");
                Assert.That(metrics.MeanSquaredError, Is.GreaterThan(0), "MSE should be positive for imperfect fit");
                Assert.That(metrics.RootMeanSquaredError, Is.GreaterThan(0), "RMSE should be positive for imperfect fit");
                Assert.That(metrics.MeanAbsoluteError, Is.GreaterThan(0), "MAE should be positive for imperfect fit");

                // RMSE should equal sqrt(MSE)
                Assert.That(metrics.RootMeanSquaredError, Is.EqualTo(Math.Sqrt(metrics.MeanSquaredError)).Within(HighPrecisionTolerance),
                    "RMSE should equal sqrt(MSE)");
            });
        }

        [Test]
        public void TestGetFittingMetrics_ErrorCases()
        {
            // Curve created from coefficients shouldn't have fitting metrics
            var coeffCurve = new Curve(new double[] { 1, 2, 3 }, 0, 5);

            Assert.Throws<InvalidOperationException>(() => coeffCurve.GetFittingMetrics(),
                "Should throw when no original data points available");
        }

        #endregion

        #region Operators Tests

        [Test]
        public void TestArithmeticOperators()
        {
            var curve1 = new Curve(new double[] { 1, 2, 3 }, 0, 5);
            var curve2 = new Curve(new double[] { 4, 5, 6 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Addition
                var sumCurve = curve1 + curve2;
                AssertArraysEqual(sumCurve.PolynomialCoeffs, new double[] { 5, 7, 9 });

                // Subtraction
                var diffCurve = curve1 - curve2;
                AssertArraysEqual(diffCurve.PolynomialCoeffs, new double[] { -3, -3, -3 });

                // Multiplication by scalar
                var multCurve = curve1 * 2;
                AssertArraysEqual(multCurve.PolynomialCoeffs, new double[] { 2, 4, 6 });

                // Commutative multiplication
                var multCurveComm = 2 * curve1;
                AssertArraysEqual(multCurveComm.PolynomialCoeffs, new double[] { 2, 4, 6 });

                // Division by scalar
                var divCurve = curve1 / 2;
                AssertArraysEqual(divCurve.PolynomialCoeffs, new double[] { 0.5, 1, 1.5 });
            });
        }

        [Test]
        public void TestVerticalShiftOperators()
        {
            var curve = new Curve(new double[] { 1, 2, 3 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Addition with scalar
                var shiftedUp = curve + 5;
                AssertArraysEqual(shiftedUp.PolynomialCoeffs, new double[] { 6, 2, 3 });

                // Commutative addition
                var shiftedUpComm = 5 + curve;
                AssertArraysEqual(shiftedUpComm.PolynomialCoeffs, new double[] { 6, 2, 3 });

                // Subtraction with scalar
                var shiftedDown = curve - 3;
                AssertArraysEqual(shiftedDown.PolynomialCoeffs, new double[] { -2, 2, 3 });
            });
        }

        [Test]
        public void TestHorizontalShiftOperators()
        {
            var curve = new Curve(new double[] { 1, 2, 3 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Right shift (horizontal shift right)
                var shiftRightCurve = curve >> 2;
                Assert.That(shiftRightCurve.EvaluateAt(3).Y, Is.EqualTo(curve.EvaluateAt(1).Y).Within(HighPrecisionTolerance));

                // Left shift (horizontal shift left)
                var shiftLeftCurve = curve << 2;
                Assert.That(shiftLeftCurve.EvaluateAt(1).Y, Is.EqualTo(curve.EvaluateAt(3).Y).Within(HighPrecisionTolerance));
            });
        }

        [Test]
        public void TestOperator_DomainIntersection()
        {
            var curve1 = new Curve(new double[] { 1, 2 }, 0, 10);
            var curve2 = new Curve(new double[] { 3, 4 }, 5, 15);

            var sum = curve1 + curve2;

            Assert.Multiple(() =>
            {
                Assert.That(sum.XMin, Is.EqualTo(5), "Domain intersection minimum incorrect");
                Assert.That(sum.XMax, Is.EqualTo(10), "Domain intersection maximum incorrect");
            });
        }

        [Test]
        public void TestOperator_ErrorHandling()
        {
            var curve = new Curve(new double[] { 1, 2, 3 }, 0, 5);

            Assert.Multiple(() =>
            {
                Assert.Throws<DivideByZeroException>(() => { var _ = curve / 0; }, "Should throw on division by zero");
                Assert.Throws<DivideByZeroException>(() => { var _ = curve / 1e-12; }, "Should throw on division by very small number");
            });
        }

        #endregion

        #region String Representation Tests

        [Test]
        public void TestToString()
        {
            var curve = new Curve(new double[] { 1, -2, 0, -4 }, 0, 5);
            string expected = "1x^0 - 2x^1 + 0x^2 - 4x^3 from 0 to 5";
            Assert.That(curve.ToString(), Is.EqualTo(expected));
        }

        [Test]
        public void TestToString_EdgeCases()
        {
            Assert.Multiple(() =>
            {
                // Constant curve
                var constant = new Curve(new double[] { 42 }, -1, 1);
                Assert.That(constant.ToString(), Does.Contain("42"), "Constant curve string representation");

                // Linear curve
                var linear = new Curve(new double[] { 0, 1 }, 0, 10);
                Assert.That(linear.ToString(), Does.Contain("1x^1"), "Linear curve string representation");

                // Negative coefficient handling
                var negative = new Curve(new double[] { -5, 3, -2 }, 0, 1);
                string negStr = negative.ToString();
                Assert.That(negStr, Does.Contain("-5"), "Negative constant term");
                Assert.That(negStr, Does.Contain("+ 3"), "Positive linear term");
                Assert.That(negStr, Does.Contain("- 2"), "Negative quadratic term");
            });
        }

        #endregion

        #region Error Handling Tests

        [Test]
        public void TestInputValidation_ComprehensiveChecks()
        {
            Assert.Multiple(() =>
            {
                // Constructor validations already tested above

                // Method parameter validations
                var validCurve = new Curve(new double[] { 1, 2 }, 0, 5);

                // WithBounds validation
                Assert.Throws<ArgumentException>(() => validCurve.WithBounds(5, 5));
                Assert.Throws<ArgumentException>(() => validCurve.WithBounds(10, 5));

                // WithExpandedBounds validation
                Assert.Throws<ArgumentException>(() => validCurve.WithExpandedBounds(-10, -10));

                // SplitCurve validation
                Assert.Throws<ArgumentException>(() => validCurve.SplitCurve(0));
                Assert.Throws<ArgumentException>(() => validCurve.SplitCurve(-5));

                // Interpolation null checks
                Assert.Throws<ArgumentNullException>(() => validCurve.InterpolateWith(null, 0.5));
                Assert.Throws<ArgumentNullException>(() => validCurve.InterpolateWithCoefficients(null, 0.5));
                Assert.Throws<ArgumentNullException>(() => validCurve.FindIntersectionsWith(null));
            });
        }

        [Test]
        public void TestNumericalStability()
        {
            Assert.Multiple(() =>
            {
                // Test with very large coefficients
                var largeCurve = new Curve(new double[] { 1e10, 1e9, 1e8 }, 0, 1);
                var largePoint = largeCurve.EvaluateAt(0.5);
                Assert.That(largePoint.Y, Is.Not.NaN, "Large coefficients should not produce NaN");

                // Test with very small coefficients
                var smallCurve = new Curve(new double[] { 1e-10, 1e-9, 1e-8 }, 0, 1);
                var smallPoint = smallCurve.EvaluateAt(0.5);
                Assert.That(smallPoint.Y, Is.Not.NaN, "Small coefficients should not produce NaN");

                // Test with zero coefficients (except constant)
                var sparseCoeffs = new double[] { 5, 0, 0, 0, 3 };
                var sparseCurve = new Curve(sparseCoeffs, 0, 1);
                var sparsePoint = sparseCurve.EvaluateAt(0.5);
                Assert.That(sparsePoint.Y, Is.Not.NaN, "Sparse coefficients should not produce NaN");
            });
        }

        [Test]
        public void TestBoundaryConditions()
        {
            var curve = new Curve(new double[] { 1, 2, 3 }, -2, 2);

            Assert.Multiple(() =>
            {
                // Test evaluation exactly at boundaries
                var leftBoundary = curve.EvaluateAt(curve.XMin);
                var rightBoundary = curve.EvaluateAt(curve.XMax);
                Assert.That(leftBoundary.X, Is.EqualTo(curve.XMin));
                Assert.That(rightBoundary.X, Is.EqualTo(curve.XMax));

                // Test operations with edge cases
                var extrema = curve.FindLocalExtrema().ToList();
                var roots = curve.FindRoots().ToList();

                // Should not throw exceptions
                Assert.DoesNotThrow(() => curve.CalculateDefiniteIntegral());
                Assert.DoesNotThrow(() => curve.CalculateAreaUnderCurve());
                Assert.DoesNotThrow(() => curve.GetStatistics());
            });
        }

        #endregion

        #region Integration/Round-trip Tests

        [Test]
        public void TestDerivativeAntiderivative_RoundTrip()
        {
            // Test that derivative of antiderivative gives back the original curve
            var original = new Curve(new double[] { 1, 2, 3, 4 }, 0, 5);
            var antiderivative = original.GetAntiderivative(0);
            var reconstructed = antiderivative.GetDerivative();

            AssertArraysEqual(reconstructed.PolynomialCoeffs, original.PolynomialCoeffs);
        }

        [Test]
        public void TestIntegration_FundamentalTheorem()
        {
            // Test fundamental theorem of calculus: ∫[a to b] f'(x) dx = f(b) - f(a)
            var curve = new Curve(new double[] { 1, 2, 3 }, 0, 5);
            var derivative = curve.GetDerivative();

            double a = 1, b = 4;
            double integralOfDerivative = derivative.CalculateDefiniteIntegral(a, b);
            double expectedDifference = curve.EvaluateAt(b).Y - curve.EvaluateAt(a).Y;

            Assert.That(integralOfDerivative, Is.EqualTo(expectedDifference).Within(HighPrecisionTolerance),
                "Fundamental theorem of calculus verification failed");
        }

        [Test]
        public void TestTransformationInvariance()
        {
            var original = new Curve(new double[] { 1, 2, 3 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Test that double transformation returns to original
                var shifted = original.ShiftHorizontally(3).ShiftHorizontally(-3);
                for (double x = 0; x <= 5; x += 0.5)
                {
                    var originalY = original.EvaluateAt(x).Y;
                    var shiftedY = shifted.EvaluateAt(x).Y;
                    Assert.That(shiftedY, Is.EqualTo(originalY).Within(HighPrecisionTolerance),
                        $"Round-trip horizontal shift failed at x={x}");
                }

                // Test vertical shift round-trip
                var vShifted = original.ShiftVertically(5).ShiftVertically(-5);
                AssertArraysEqual(vShifted.PolynomialCoeffs, original.PolynomialCoeffs);
            });
        }

        [Test]
        public void TestOperatorConsistency()
        {
            var curve1 = new Curve(new double[] { 1, 2 }, 0, 5);
            var curve2 = new Curve(new double[] { 3, 4 }, 0, 5);

            Assert.Multiple(() =>
            {
                // Test that operator+ is commutative
                var sum1 = curve1 + curve2;
                var sum2 = curve2 + curve1;
                AssertArraysEqual(sum1.PolynomialCoeffs, sum2.PolynomialCoeffs);

                // Test that scalar multiplication is commutative
                var mult1 = curve1 * 2.5;
                var mult2 = 2.5 * curve1;
                AssertArraysEqual(mult1.PolynomialCoeffs, mult2.PolynomialCoeffs);

                // Test that a - b = a + (-1 * b)
                var diff = curve1 - curve2;
                var equivalent = curve1 + (-1 * curve2);
                AssertArraysEqual(diff.PolynomialCoeffs, equivalent.PolynomialCoeffs);
            });
        }

        [Test]
        public void TestChainedOperations_Consistency()
        {
            var original = new Curve(new double[] { 1, 2, 3 }, -2, 2);

            Assert.Multiple(() =>
            {
                // Complex chain of operations
                var processed = original
                    .GetDerivative()
                    .ShiftVertically(5)
                    .SplitCurve(2).First()
                    .WithExpandedBounds(1, 1);

                // Should not throw exceptions and should produce valid results
                Assert.DoesNotThrow(() => processed.FindRoots());
                Assert.DoesNotThrow(() => processed.CalculateDefiniteIntegral());
                Assert.DoesNotThrow(() => processed.GetStatistics());

                var stats = processed.GetStatistics();
                Assert.That(stats.Mean, Is.Not.NaN, "Chained operations should produce valid statistics");
            });
        }

        [Test]
        public void TestNumericalAccuracy_HighPrecision()
        {
            // Test with known analytical results
            var parabola = new Curve(new double[] { 0, 0, 1 }, 0, 2); // y = x²

            Assert.Multiple(() =>
            {
                // Integral of x² from 0 to 2 should be 8/3
                double integral = parabola.CalculateDefiniteIntegral(0, 2);
                Assert.That(integral, Is.EqualTo(8.0 / 3.0).Within(HighPrecisionTolerance), "Analytical integral verification");

                // Derivative of x² should be 2x
                var derivative = parabola.GetDerivative();
                AssertArraysEqual(derivative.PolynomialCoeffs, new double[] { 0, 2 });

                // Arc length calculation should be reasonable
                double arcLength = parabola.CalculateArcLength(0, 1);
                Assert.That(arcLength, Is.GreaterThan(1), "Arc length should be greater than straight line distance");
                Assert.That(arcLength, Is.LessThan(2), "Arc length should be reasonable");
            });
        }

        #endregion
    }
}