from math import isclose

from material.integrate import QuadpackIntegrator

from .test_IntegratorInterface import ParabolicTestIntegrand

def test_quadpack_integrator():
    integrator = QuadpackIntegrator()

    parabolic_test = ParabolicTestIntegrand()
    parabolic_result = integrator.integrate(parabolic_test)
    assert isclose(parabolic_result, parabolic_test.analytical())
