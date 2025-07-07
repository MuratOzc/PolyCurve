using PolyCurve;

namespace Tests
{
    public abstract class BaseCurveTests
    {
        protected const double Tolerance = 0.15; // Some methods can handle much smalletr tolerances, but this is a good default for most tests

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

        [Test]
        public void TestArcLength()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateArcLength(-3, 11), Is.EqualTo(2598.79744060024886176663).Within(Tolerance), "Polynomial curve arc length failed");
                Assert.That(sinusoidalCurve.CalculateArcLength(0, 2 * Math.PI), Is.EqualTo(21.393860164471).Within(Tolerance), "Sinusoidal curve arc length failed");
                Assert.That(exponentialCurve.CalculateArcLength(0, 5), Is.EqualTo(31.8893).Within(Tolerance), "Exponential curve arc length failed");
            });
        }

        [Test]
        public void TestArea()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateAreaUnderCurve(-3, 11), Is.EqualTo(2561.2).Within(Tolerance), "Polynomial curve area failed");
                Assert.That(sinusoidalCurve.CalculateAreaUnderCurve(0, 2 * Math.PI), Is.EqualTo(21.6224).Within(Tolerance), "Sinusoidal curve area failed");
                Assert.That(exponentialCurve.CalculateAreaUnderCurve(0, 5), Is.EqualTo(44.724).Within(Tolerance), "Exponential curve area failed");
            });
        }

        [Test]
        public void TestIntegral()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateDefiniteIntegral(-3, 11), Is.EqualTo(1178.8).Within(Tolerance), "Polynomial curve integral failed");
                Assert.That(sinusoidalCurve.CalculateDefiniteIntegral(0, 2 * Math.PI), Is.EqualTo(12.5664).Within(Tolerance), "Sinusoidal curve integral failed");
                Assert.That(exponentialCurve.CalculateDefiniteIntegral(0, 5), Is.EqualTo(44.724).Within(Tolerance), "Exponential curve integral failed");
            });
        }

        [Test]
        public void TestCurvature()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateCurvatureAt(0), Is.EqualTo(0.000457596).Within(Tolerance), "Polynomial curve curvature failed");
                Assert.That(sinusoidalCurve.CalculateCurvatureAt(Math.PI / 2), Is.EqualTo(5).Within(Tolerance), "Sinusoidal curve curvature failed");
                Assert.That(exponentialCurve.CalculateCurvatureAt(1), Is.EqualTo(0.1924).Within(Tolerance), "Exponential curve curvature failed");
            });
        }

        [Test]
        public void TestInflectionPoints()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindInflectionPoints(), [new Point2D(1.55051, -66.3837), new Point2D(6.44949, 90.3837)]);
                AssertPointsEqual(sinusoidalCurve.FindInflectionPoints(), [new Point2D(0, 2.0), new Point2D(Math.PI, 2.0), new Point2D(2 * Math.PI, 2.0)]);
                Assert.That(exponentialCurve.FindInflectionPoints(), Is.Empty, "Exponential curve inflection points failed");

            });
        }

        [Test]
        public void TestLocalExtrema()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindLocalExtrema(), new[] { new Point2D(-0.44949, -271.15), new Point2D(4.4495, 199.15), new Point2D(8, 0) });
                AssertPointsEqual(sinusoidalCurve.FindLocalExtrema(), new[] { new Point2D(Math.PI / 2, 7.0), new Point2D(3 * Math.PI / 2, -3) });
                Assert.That(exponentialCurve.FindLocalExtrema(), Is.Empty, "Exponential curve local extrema failed");
            });
        }

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
                //AssertPointEqual(polynomialCurve.FindMinimum(), (-0.44949, -271.15));
                //AssertPointEqual(sinusoidalCurve.FindMinimum(), (3 * Math.PI / 2, -3));
                //AssertPointEqual(exponentialCurve.FindMinimum(), (0, 1));

                AssertPointEqual(polynomialCurve.FindMinimum(), new Point2D(-0.44949, -271.15));
                AssertPointEqual(sinusoidalCurve.FindMinimum(), new Point2D(3 * Math.PI / 2, -3));
                AssertPointEqual(exponentialCurve.FindMinimum(), new Point2D(0, 1));

            });
        }

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
        public void TestFindXForY()
        {
            Assert.Multiple(() =>
            {
                AssertArraysEqual(polynomialCurve.FindXValuesForY(100), new[] { -2.2263, 2.7659, 6.3370, 9.1234 });
                AssertArraysEqual(sinusoidalCurve.FindXValuesForY(4), new[] { 0.4636, 2.6780 });
                AssertArraysEqual(exponentialCurve.FindXValuesForY(8), new[] { 3.0 });
            });
        }

        [Test]
        public void TestInterpolate()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.EvaluateAt(5).Y, Is.EqualTo(189.0).Within(Tolerance), "Polynomial curve interpolation failed");
                Assert.That(sinusoidalCurve.EvaluateAt(Math.PI / 4).Y, Is.EqualTo(5.5355).Within(Tolerance), "Sinusoidal curve interpolation failed");
                Assert.That(exponentialCurve.EvaluateAt(2).Y, Is.EqualTo(4.0).Within(Tolerance), "Exponential curve interpolation failed");
            });
        }

        [Test]
        public void TestNormalLine()
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
        public void TestTangentLine()
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
        public void TestDerivativeCurve()
        {
            Assert.Multiple(() =>
            {
                var polyDerivative = polynomialCurve.GetDerivative();
                AssertArraysEqual(polyDerivative.PolynomialCoeffs, new[] { 64.0, 120, -48, 4 });

                var sinDerivative = sinusoidalCurve.GetDerivative();
                Assert.That(sinDerivative.EvaluateAt(Math.PI / 2).Y, Is.EqualTo(0).Within(Tolerance), "Sinusoidal derivative at PI/2 failed");
                Assert.That(sinDerivative.EvaluateAt(0).Y, Is.EqualTo(5).Within(Tolerance), "Sinusoidal derivative at 0 failed");

                var expDerivative = exponentialCurve.GetDerivative();
                Assert.That(expDerivative.EvaluateAt(1).Y, Is.EqualTo(1.3863).Within(Tolerance), "Exponential derivative failed");
            });
        }

        [Test]
        public void TestIntegralCurve()
        {
            Assert.Multiple(() =>
            {
                var polyIntegral = polynomialCurve.GetAntiderivative();
                Assert.That(polyIntegral.EvaluateAt(5).Y, Is.EqualTo(145).Within(Tolerance), "Polynomial integral curve failed");

                var sinIntegral = sinusoidalCurve.GetAntiderivative(integrationConstant: -5);
                Assert.That(sinIntegral.EvaluateAt(Math.PI).Y, Is.EqualTo(5 + 2 * Math.PI).Within(Tolerance), "Sinusoidal integral curve failed");

                var expIntegral = exponentialCurve.GetAntiderivative(integrationConstant: 1.4427);
                Assert.That(expIntegral.EvaluateAt(1).Y, Is.EqualTo(2.88539).Within(Tolerance), "Exponential integral curve failed");
            });
        }

        [Test]
        public void TestSplitCurve()
        {
            Assert.Multiple(() =>
            {
                Curve[] polyCurves = polynomialCurve.SplitCurve([4]).ToArray();
                Assert.That(polyCurves[0].XMax, Is.EqualTo(4).Within(Tolerance), "Polynomial left split failed");
                Assert.That(polyCurves[1].XMin, Is.EqualTo(4).Within(Tolerance), "Polynomial right split failed");

                Curve[] sinCurves = sinusoidalCurve.SplitCurve([Math.PI]).ToArray();
                Assert.That(sinCurves[0].XMax, Is.EqualTo(Math.PI).Within(Tolerance), "Sinusoidal left split failed");
                Assert.That(sinCurves[1].XMin, Is.EqualTo(Math.PI).Within(Tolerance), "Sinusoidal right split failed");

                var expCurves = exponentialCurve.SplitCurve([2.5]).ToArray();
                Assert.That(expCurves[0].XMax, Is.EqualTo(2.5).Within(Tolerance), "Exponential left split failed");
                Assert.That(expCurves[1].XMin, Is.EqualTo(2.5).Within(Tolerance), "Exponential right split failed");
            });
        }

        [Test]
        public void TestOperators()
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

                // Division by scalar
                var divCurve = curve1 / 2;
                AssertArraysEqual(divCurve.PolynomialCoeffs, new double[] { 0.5, 1, 1.5 });

                // Right shift (horizontal shift right)
                var shiftRightCurve = curve1 >> 2;
                Assert.That(shiftRightCurve.EvaluateAt(3).Y, Is.EqualTo(curve1.EvaluateAt(1).Y).Within(Tolerance));

                // Left shift (horizontal shift left)
                var shiftLeftCurve = curve1 << 2;
                Assert.That(shiftLeftCurve.EvaluateAt(1).Y, Is.EqualTo(curve1.EvaluateAt(3).Y).Within(Tolerance));
            });
        }

        [Test]
        public void TestFindIntersections()
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
        public void TestToString()
        {
            var curve = new Curve(new double[] { 1, -2, 0, -4 }, 0, 5);
            string expected = "1x^0 - 2x^1 + 0x^2 - 4x^3 from 0 to 5";
            Assert.That(curve.ToString(), Is.EqualTo(expected));
        }

        // NEW DOMAIN OPERATION TESTS BELOW

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

        [Test]
        public void TestDomainOperations_ExtremaAnalysis()
        {
            // Test curve y = x^2 - 4 has minimum at (0, -4)
            var originalExtrema = testCurve.FindLocalExtrema().ToList();
            var originalMin = testCurve.FindMinimum();

            // Contract domain to exclude the minimum
            var contractedExtrema = testCurve.WithBounds(1, 3).FindLocalExtrema().ToList();
            var contractedMin = testCurve.WithBounds(1, 3).FindMinimum();

            Assert.Multiple(() =>
            {
                Assert.That(originalExtrema.Count, Is.EqualTo(1), "Should find 1 extremum in original domain");
                AssertPointEqual(originalMin, new Point2D(0, -4));

                Assert.That(contractedExtrema.Count, Is.EqualTo(0), "Should find no local extrema in contracted domain");
                AssertPointEqual(contractedMin, new Point2D(1, -3)); // Minimum at left boundary
            });
        }

        [Test]
        public void TestDomainOperations_Integration()
        {
            // Test that integration works correctly with different bounds
            double originalIntegral = testCurve.CalculateDefiniteIntegral();
            double expandedIntegral = testCurve.WithExpandedBounds(1, 1).CalculateDefiniteIntegral();
            double contractedIntegral = testCurve.WithBounds(-1, 1).CalculateDefiniteIntegral();

            // For y = x^2 - 4, integral from -1 to 1 = [x^3/3 - 4x] = (1/3 - 4) - (-1/3 + 4) = -22/3
            double expectedContracted = -22.0 / 3.0;

            Assert.Multiple(() =>
            {
                Assert.That(originalIntegral, Is.Not.EqualTo(expandedIntegral), "Expanded integral should differ");
                Assert.That(contractedIntegral, Is.EqualTo(expectedContracted).Within(Tolerance), "Contracted integral incorrect");
            });
        }

        [Test]
        public void TestDomainOperations_ChainedOperations()
        {
            // Test chaining operations with temporary bounds
            var chainedResults = testCurve
                .WithExpandedBounds(2, 2)
                .FindLocalExtrema()
                .Where(p => p.Y < 0)
                .ToList();

            var directResults = testCurve.FindLocalExtrema().Where(p => p.Y < 0).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(chainedResults.Count, Is.EqualTo(directResults.Count), "Chained results should match direct results");
                if (chainedResults.Any() && directResults.Any())
                {
                    AssertPointEqual(chainedResults.First(), directResults.First());
                }
            });
        }

        [Test]
        public void TestDomainOperations_Statistics()
        {
            var originalStats = testCurve.GetStatistics();
            var expandedStats = testCurve.WithExpandedBounds(1, 1).GetStatistics();
            var contractedStats = testCurve.WithBounds(-1, 1).GetStatistics();

            Assert.Multiple(() =>
            {
                Assert.That(originalStats.Mean, Is.Not.EqualTo(expandedStats.Mean), "Expanded stats should differ");
                Assert.That(originalStats.Mean, Is.Not.EqualTo(contractedStats.Mean), "Contracted stats should differ");
                Assert.That(expandedStats.CentroidPoint.X, Is.EqualTo(0).Within(Tolerance), "Centroid should be at x=0 for symmetric curve");
            });
        }

        [Test]
        public void TestDomainOperations_Intersections()
        {
            // Create a line y = 0 (horizontal line through roots)
            var horizontalLine = new Curve(new double[] { 0 }, -10, 10);

            // Find intersections with different domain bounds
            var originalIntersections = testCurve.FindIntersectionsWith(horizontalLine).ToList();
            var expandedIntersections = testCurve.WithBounds(-5, 5)
                .FindIntersectionsWith(horizontalLine.WithBounds(-5, 5)).ToList();
            var contractedIntersections = testCurve.WithBounds(-1, 1)
                .FindIntersectionsWith(horizontalLine.WithBounds(-1, 1)).ToList();

            Assert.Multiple(() =>
            {
                Assert.That(originalIntersections.Count, Is.EqualTo(2), "Should find 2 intersections in original domain");
                Assert.That(expandedIntersections.Count, Is.EqualTo(2), "Should find 2 intersections in expanded domain");
                Assert.That(contractedIntersections.Count, Is.EqualTo(0), "Should find no intersections in contracted domain");
            });
        }
    }
}