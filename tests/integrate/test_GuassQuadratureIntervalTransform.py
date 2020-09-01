from math import isclose

from super_material.integrate import (
    QuadpackIntegrator,
    LinearGuassQuadratureIntervalTransform,
    TransformedIntegrand,
)

from .test_IntegratorInterface import ParabolicTestIntegrand


def test_linear_guass_quadrature_integrand_transform():
    integrator = QuadpackIntegrator()

    # Parabolic test
    parabolic_test = ParabolicTestIntegrand()
    parabolic_transform = LinearGuassQuadratureIntervalTransform(
        parabolic_test.interval()
    )
    transformed_test = TransformedIntegrand(parabolic_test, parabolic_transform)

    # Interval
    transformed_interval = transformed_test.interval()
    assert transformed_interval.start().defined_on_boundary()
    assert transformed_interval.end().defined_on_boundary()
    assert isclose(transformed_interval.start().value(), -1)
    assert isclose(transformed_interval.end().value(), 1)

    # Integration
    transformed_result = integrator.integrate(transformed_test)
    assert isclose(transformed_result, parabolic_test.analytical())
