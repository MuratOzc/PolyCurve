using CurvesLibrary;
using MathNet.Numerics;

namespace Tests
{
    public abstract class BaseCurveTests
    {
        protected Curve curve;
        protected const double Tolerance = 0.15;

        protected void AssertPointsEqual(IEnumerable<(double X, double Y)> actual, IEnumerable<(double X, double Y)> expected)
        {
            var actualList = actual.OrderBy(p => p.X).ToList();
            var expectedList = expected.OrderBy(p => p.X).ToList();
            Assert.That(actualList.Count, Is.EqualTo(expectedList.Count));
            for (int i = 0; i < actualList.Count; i++)
            {
                AssertPointEqual(actualList[i], expectedList[i]);
            }
        }

        protected void AssertPointEqual((double X, double Y) actual, (double X, double Y) expected)
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
    }

    [TestFixture]
    public class CurveTests : BaseCurveTests
    {
        private Curve polynomialCurve;
        private Curve sinusoidalCurve;
        private Curve exponentialCurve;

        [SetUp]
        public void Setup()
        {
            // Polynomial curve: y = -256 + 64x + 60x^2 - 16x^3 + x^4
            double[] coeffs = { -256, 64, 60, -16, 1 };
            polynomialCurve = new Curve(coeffs, -3, 11);

            // Sinusoidal curve: y = 5 * sin(x) + 2
            Func<double, double> sinFunc = x => 5 * Math.Sin(x) + 2;
            double[] xValues = Generate.LinearSpaced(100, 0, 2 * Math.PI);
            double[] yValues = xValues.Select(sinFunc).ToArray();
            sinusoidalCurve = new Curve(xValues, yValues, polynomialDegree: 10);

            // Exponential curve: y = 2^x
            Func<double, double> expFunc = x => Math.Pow(2, x);
            xValues = Generate.LinearSpaced(100, 0, 5);
            yValues = xValues.Select(expFunc).ToArray();
            exponentialCurve = new Curve(xValues, yValues, polynomialDegree: 10);
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
                Assert.That(polynomialCurve.CalculateArea(-3, 11), Is.EqualTo(2561.2).Within(Tolerance), "Polynomial curve area failed");
                Assert.That(sinusoidalCurve.CalculateArea(0, 2 * Math.PI), Is.EqualTo(21.6224).Within(Tolerance), "Sinusoidal curve area failed");
                Assert.That(exponentialCurve.CalculateArea(0, 5), Is.EqualTo(44.724).Within(Tolerance), "Exponential curve area failed");
            });
        }

        [Test]
        public void TestIntegral()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateIntegral(-3, 11), Is.EqualTo(1178.8).Within(Tolerance), "Polynomial curve integral failed");
                Assert.That(sinusoidalCurve.CalculateIntegral(0, 2 * Math.PI), Is.EqualTo(12.5664).Within(Tolerance), "Sinusoidal curve integral failed");
                Assert.That(exponentialCurve.CalculateIntegral(0, 5), Is.EqualTo(44.724).Within(Tolerance), "Exponential curve integral failed");
            });
        }

        [Test]
        public void TestCurvature()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.CalculateCurvature(0), Is.EqualTo(0.000457596).Within(Tolerance), "Polynomial curve curvature failed");
                Assert.That(sinusoidalCurve.CalculateCurvature(Math.PI / 2), Is.EqualTo(5).Within(Tolerance), "Sinusoidal curve curvature failed");
                Assert.That(exponentialCurve.CalculateCurvature(1), Is.EqualTo(0.1924).Within(Tolerance), "Exponential curve curvature failed");
            });
        }

        [Test]
        public void TestInflectionPoints()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindInflectionPoints(), new[] { (1.55051, -66.3837), (6.44949, 90.3837) });
                AssertPointsEqual(sinusoidalCurve.FindInflectionPoints(), new[] { (0, 2.0), (Math.PI, 2.0), (2 * Math.PI, 2.0) });
                Assert.That(exponentialCurve.FindInflectionPoints(), Is.Empty, "Exponential curve inflection points failed");
            });
        }

        [Test]
        public void TestLocalExtrema()
        {
            Assert.Multiple(() =>
            {
                AssertPointsEqual(polynomialCurve.FindLocalExtrema(), new[] { (-0.44949, -271.15), (4.4495, 199.15), (8, 0) });
                AssertPointsEqual(sinusoidalCurve.FindLocalExtrema(), new[] { (Math.PI / 2, 7.0), (3 * Math.PI / 2, -3) });
                Assert.That(exponentialCurve.FindLocalExtrema(), Is.Empty, "Exponential curve local extrema failed");
            });
        }

        [Test]
        public void TestFindMaximum()
        {
            Assert.Multiple(() =>
            {
                AssertPointEqual(polynomialCurve.FindMaximum(), (11.0, 1053.0));
                AssertPointEqual(sinusoidalCurve.FindMaximum(), (Math.PI / 2, 7));
                AssertPointEqual(exponentialCurve.FindMaximum(), (5, 32));
            });
        }

        [Test]
        public void TestFindMinimum()
        {
            Assert.Multiple(() =>
            {
                AssertPointEqual(polynomialCurve.FindMinimum(), (-0.44949, -271.15));
                AssertPointEqual(sinusoidalCurve.FindMinimum(), (3 * Math.PI / 2, -3));
                AssertPointEqual(exponentialCurve.FindMinimum(), (0, 1));
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
                AssertArraysEqual(polynomialCurve.FindXForY(100), new[] { -2.2263, 2.7659, 6.3370, 9.1234 });
                AssertArraysEqual(sinusoidalCurve.FindXForY(4), new[] { 0.4636, 2.6780 });
                AssertArraysEqual(exponentialCurve.FindXForY(8), new[] { 3.0 });
            });
        }

        [Test]
        public void TestInterpolate()
        {
            Assert.Multiple(() =>
            {
                Assert.That(polynomialCurve.Interpolate(5), Is.EqualTo(189.0).Within(Tolerance), "Polynomial curve interpolation failed");
                Assert.That(sinusoidalCurve.Interpolate(Math.PI / 4), Is.EqualTo(5.5355).Within(Tolerance), "Sinusoidal curve interpolation failed");
                Assert.That(exponentialCurve.Interpolate(2), Is.EqualTo(4.0).Within(Tolerance), "Exponential curve interpolation failed");
            });
        }

        [Test]
        public void TestNormalLine()
        {
            Assert.Multiple(() =>
            {
                var polyNormal = polynomialCurve.NormalLine(5);
                AssertArraysEqual(polyNormal.PolynomialCoeffs, new[] { 188.861, 0.0277778 });

                var sinNormal = sinusoidalCurve.NormalLine(Math.PI / 4);
                AssertArraysEqual(sinNormal.PolynomialCoeffs, new[] { 5.75768, -0.282843 });

                var expNormal = exponentialCurve.NormalLine(2);
                AssertArraysEqual(expNormal.PolynomialCoeffs, new[] { 4.72135, -0.360674 });
            });
        }

        [Test]
        public void TestTangentLine()
        {
            Assert.Multiple(() =>
            {
                var polyTangent = polynomialCurve.TangentLine(5);
                AssertArraysEqual(polyTangent.PolynomialCoeffs, new[] { 369.0, -36.0 });

                var sinTangent = sinusoidalCurve.TangentLine(Math.PI / 4);
                AssertArraysEqual(sinTangent.PolynomialCoeffs, new[] { 2.75873, 3.53553 });

                var expTangent = exponentialCurve.TangentLine(2);
                AssertArraysEqual(expTangent.PolynomialCoeffs, new[] { -1.54518, 2.77259 });
            });
        }

        [Test]
        public void TestDerivativeCurve()
        {
            Assert.Multiple(() =>
            {
                var polyDerivative = polynomialCurve.DerivativeCurve();
                AssertArraysEqual(polyDerivative.PolynomialCoeffs, new[] { 64.0, 120, -48, 4 });

                var sinDerivative = sinusoidalCurve.DerivativeCurve();
                Assert.That(sinDerivative.Interpolate(Math.PI / 2), Is.EqualTo(0).Within(Tolerance), "Sinusoidal derivative at PI/2 failed");
                Assert.That(sinDerivative.Interpolate(0), Is.EqualTo(5).Within(Tolerance), "Sinusoidal derivative at 0 failed");

                var expDerivative = exponentialCurve.DerivativeCurve();
                Assert.That(expDerivative.Interpolate(1), Is.EqualTo(1.3863).Within(Tolerance), "Exponential derivative failed");
            });
        }

        [Test]
        public void TestIntegralCurve()
        {
            Assert.Multiple(() =>
            {
                var polyIntegral = polynomialCurve.IntegralCurve();
                Assert.That(polyIntegral.Interpolate(5), Is.EqualTo(145).Within(Tolerance), "Polynomial integral curve failed");

                var sinIntegral = sinusoidalCurve.IntegralCurve(c: -5);
                Assert.That(sinIntegral.Interpolate(Math.PI), Is.EqualTo(5 + 2 * Math.PI).Within(Tolerance), "Sinusoidal integral curve failed");

                var expIntegral = exponentialCurve.IntegralCurve(c: 1.4427);
                Assert.That(expIntegral.Interpolate(1), Is.EqualTo(2.88539).Within(Tolerance), "Exponential integral curve failed");
            });
        }

        [Test]
        public void TestSplitCurve()
        {
            Assert.Multiple(() =>
            {
                var (polyLeft, polyRight) = polynomialCurve.SplitCurve(4);
                Assert.That(polyLeft.XMax, Is.EqualTo(4).Within(Tolerance), "Polynomial left split failed");
                Assert.That(polyRight.XMin, Is.EqualTo(4).Within(Tolerance), "Polynomial right split failed");

                var (sinLeft, sinRight) = sinusoidalCurve.SplitCurve(Math.PI);
                Assert.That(sinLeft.XMax, Is.EqualTo(Math.PI).Within(Tolerance), "Sinusoidal left split failed");
                Assert.That(sinRight.XMin, Is.EqualTo(Math.PI).Within(Tolerance), "Sinusoidal right split failed");

                var (expLeft, expRight) = exponentialCurve.SplitCurve(2.5);
                Assert.That(expLeft.XMax, Is.EqualTo(2.5).Within(Tolerance), "Exponential left split failed");
                Assert.That(expRight.XMin, Is.EqualTo(2.5).Within(Tolerance), "Exponential right split failed");
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
                Assert.That(shiftRightCurve.Interpolate(3), Is.EqualTo(curve1.Interpolate(1)).Within(Tolerance));

                // Left shift (horizontal shift left)
                var shiftLeftCurve = curve1 << 2;
                Assert.That(shiftLeftCurve.Interpolate(1), Is.EqualTo(curve1.Interpolate(3)).Within(Tolerance));
            });
        }

        [Test]
        public void TestFindIntersections()
        {
            var curve1 = new Curve(new double[] { 0, 1 }, 0, 10); // y = x
            var curve2 = new Curve(new double[] { 2, -1 }, 0, 10); // y = 2 - x

            var intersections = curve1.FindIntersections(curve2);
            var expected = new List<(double X, double Y)>
            {
                (1, 1),
            };

            var sinExpoIntersections = sinusoidalCurve.FindIntersections(exponentialCurve);
            var sinExpoExpected = new List<(double X, double Y)>
            {
                (2.41344, 5.32743),
            };
            var sinPolyIntersections = sinusoidalCurve.FindIntersections(polynomialCurve);
            var sinPolyExpected = new List<(double X, double Y)>
            {
                (2.04495, 6.44783),
            };
            var expoPolyIntersections = exponentialCurve.FindIntersections(polynomialCurve);
            var expoPolyExpected = new List<(double X, double Y)>
            {
                (2.0284, 4.07952),
            };

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
            var curve = new Curve(new double[] { 1, -2, 3, -4 }, 0, 5);
            string expected = "1x^0 + -2x^1 + 3x^2 + -4x^3";
            Assert.That(curve.ToString(), Is.EqualTo(expected));
        }
    }
}