from math import isclose, atanh, atan, sqrt, pi

from material.conductivity.MattisBardeenComplexConductivity import *


from material.integrate import (
    QuadpackIntegrator,
    TransformedIntegrand,
)


class MattisBardeenLowerSingularityTestIntegrand(IntegrandInterface):
    # \int_{a}^{b} \frac{1}{\sqrt{x^2-a^2}}
    _a: float
    _b: float

    def __init__(self, a, b):
        self._a = a
        self._b = b

    def evaluate(self, x: float) -> float:
        out = 1 / sqrt(x ** 2 - self._a ** 2)
        return out

    def interval(self) -> IntegrandInterval:
        start = IntegrandBoundary(self._a, False)
        end = IntegrandBoundary(self._b, True)
        interval = IntegrandInterval(start, end)
        return interval

    def analytical(self) -> float:
        return 2 * atanh(sqrt((self._b - self._a) / (self._a + self._b)))


class MattisBardeenUpperSingularityTestIntegrand(IntegrandInterface):
    # \int_{a}^{b} \frac{1}{\sqrt{b^2-x**2}}
    _a: float
    _b: float

    def __init__(self, a, b):
        self._a = a
        self._b = b

    def evaluate(self, x: float) -> float:
        out = 1 / sqrt(self._b ** 2 - x ** 2)
        return out

    def interval(self) -> IntegrandInterval:
        start = IntegrandBoundary(self._a, True)
        end = IntegrandBoundary(self._b, False)
        interval = IntegrandInterval(start, end)
        return interval

    def analytical(self) -> float:
        return 2 * atan(sqrt((self._b - self._a) / (self._a + self._b)))


class MattisBardeenUpperAndLowerSingularityTestIntegrand(IntegrandInterface):
    # \int_{a}^{B} \frac{1}{\sqrt{1- x^2}}

    def evaluate(self, x: float) -> float:
        return 1 / sqrt(1 - x ** 2)

    def interval(self) -> IntegrandInterval:
        start = IntegrandBoundary(-1, False)
        end = IntegrandBoundary(1, False)
        interval = IntegrandInterval(start, end)
        return interval

    def analytical(self) -> float:
        return pi


def test_mattis_bardeen_remove_lower_singularity_transform():
    a = 1
    b = 2

    integrator = QuadpackIntegrator()
    transform = MattisBardeenRemoveLowerSingularityTransform(a)

    # Integratable singularity
    lower_singular_test = MattisBardeenLowerSingularityTestIntegrand(a, b)
    transformed_test = TransformedIntegrand(lower_singular_test, transform)

    # Interval
    transformed_interval = transformed_test.interval()
    assert isclose(transformed_interval.start().value(), 0)
    assert isclose(transformed_interval.end().value(), sqrt(b - a))

    # Transformed evaluations
    assert isclose(transformed_test.evaluate(1), 2 / sqrt(3))

    # Integration
    transformed_result = integrator.integrate(transformed_test)
    assert isclose(transformed_result, lower_singular_test.analytical())


def test_mattis_bardeen_remove_upper_singularity_transform():
    a = 1
    b = 2

    integrator = QuadpackIntegrator()
    transform = MattisBardeenRemoveUpperSingularityTransform(b)

    # Integratable singularity
    upper_singular_test = MattisBardeenUpperSingularityTestIntegrand(a, b)
    transformed_test = TransformedIntegrand(upper_singular_test, transform)

    # Interval
    transformed_interval = transformed_test.interval()
    assert isclose(transformed_interval.start().value(), sqrt(b - 1))
    assert isclose(transformed_interval.end().value(), 0)

    # Transformed evaluations
    assert isclose(transformed_test.evaluate(1), -2 / sqrt(3))

    # Integration
    transformed_result = integrator.integrate(transformed_test)
    assert isclose(transformed_result, upper_singular_test.analytical())


def test_mattis_bardeen_remove_upper_and_lower_singularity_transform():
    relative_tolerance = 1e-12
    absolute_tolerance = 1e-12

    integrator = QuadpackIntegrator(
        relative_tolerance=relative_tolerance,
        absolute_tolerance=absolute_tolerance,
        limit=10,
    )
    lower_transform = MattisBardeenRemoveLowerSingularityTransform(-1)
    upper_transform = MattisBardeenRemoveUpperSingularityTransform(sqrt(2))

    integrand_test = MattisBardeenUpperAndLowerSingularityTestIntegrand()
    first_transformed = TransformedIntegrand(integrand_test, lower_transform)
    second_transformed = TransformedIntegrand(first_transformed, upper_transform)

    result = integrator.integrate(second_transformed)
    assert isclose(
        result,
        integrand_test.analytical(),
        rel_tol=relative_tolerance,
        abs_tol=absolute_tolerance,
    )
