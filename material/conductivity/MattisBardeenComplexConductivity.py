from math import sqrt, exp, pi, inf

from .SuperconductorConductivityInterface import SuperconductorConductivityInterface

from ..gap_energy.GapEnergyInterface import GapEnergyInterface
from ..integrate import (
    IntegrandBoundary,
    IntegrandInterface,
    IntegrandInterval,
    IntegrandIntervalTransformInterface,
    QuadpackIntegrator,
    TransformedIntegrand,
)

from ..constants import h_bar, k_B


def f(E, T):
    exponent = E / (k_B * T)

    if exponent > 500:
        return 0

    denumerator = exp(exponent) + 1
    return 1 / denumerator


class MattisBardeenRemoveLowerSingularityTransform(IntegrandIntervalTransformInterface):
    _singularity: float

    def __init__(self, singularity: float):
        self._singularity = singularity

    def transform(self, x: float) -> float:
        return sqrt(x - self._singularity)

    def inverse_transform(self, u: float) -> float:
        return self._singularity + u ** 2

    def transform_derivative(self, x: float) -> float:
        return 0.5 / sqrt(x - self._singularity)


class MattisBardeenRemoveUpperSingularityTransform(IntegrandIntervalTransformInterface):
    _singularity: float

    def __init__(self, singularity: float):
        self._singularity = singularity

    def transform(self, x: float) -> float:
        return sqrt(self._singularity - x)

    def inverse_transform(self, u: float) -> float:
        return self._singularity - u ** 2

    def transform_derivative(self, x: float) -> float:
        return -0.5 / sqrt(self._singularity - x)


def remove_lower_singularity(integrand: IntegrandInterface):
    lower_singularity = integrand.interval().start().value()
    transform = MattisBardeenRemoveLowerSingularityTransform(lower_singularity)
    transformed_integrand = TransformedIntegrand(integrand, transform)
    return transformed_integrand


def remove_upper_singularity(integrand: IntegrandInterface):
    upper_singularity = integrand.interval().end().value()
    transform = MattisBardeenRemoveUpperSingularityTransform(upper_singularity)
    transformed_integrand = TransformedIntegrand(integrand, transform)
    return transformed_integrand


class MattisBardeenRealFirstIntegrand(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float

    def __init__(self, gap_energy: float, temperature: float, omega: float):
        self._gap_energy = gap_energy
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(self._gap_energy, False)
        upper = IntegrandBoundary(inf, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        a = f(E, self._temperature) - f(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        c = sqrt(E ** 2 - self._gap_energy ** 2)
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy ** 2)

        return a * b / (c * d)


class MattisBardeenRealSecondIntegrand(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float

    def __init__(self, gap_energy: float, temperature: float, omega: float):
        assert h_bar * omega > 2 * gap_energy

        self._gap_energy = gap_energy
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(self._gap_energy - h_bar * self._omega, False)
        upper = IntegrandBoundary(-self._gap_energy, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        a = 1 - 2 * f(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        c = sqrt(E ** 2 - self._gap_energy ** 2)
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy ** 2)

        return a * b / (c * d)


class MattisBardeenImaginaryIntegrand(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float

    def __init__(self, gap_energy: float, temperature: float, omega: float):
        self._gap_energy = gap_energy
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = self._gap_energy - h_bar * self._omega

        if lower > -self._gap_energy:
            start = IntegrandBoundary(lower, False)
        else:
            start = IntegrandBoundary(-self._gap_energy, True)

        end = IntegrandBoundary(self._gap_energy, False)

        return IntegrandInterval(start, end)

    def evaluate(self, E: float) -> float:
        a = 1 - 2 * f(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        c = sqrt(self._gap_energy ** 2 - E ** 2)
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy ** 2)

        return a * b / (c * d)


class MattisBardeenComplexConductivity(SuperconductorConductivityInterface):
    _gap_energy: GapEnergyInterface
    _conductivity_0: float
    _integrator: QuadpackIntegrator

    def __init__(self, gap_energy: GapEnergyInterface, conductivity_0: float):
        self._gap_energy = gap_energy
        self._conductivity_0 = conductivity_0
        self._integrator = QuadpackIntegrator(
            absolute_tolerance=1e-12, relative_tolerance=1e-12, limit=10
        )

    def evaluate_first_real_integral(
        self, gap_energy: float, temperature: float, omega: float
    ) -> float:
        integrand = MattisBardeenRealFirstIntegrand(gap_energy, temperature, omega)
        integrand = remove_lower_singularity(integrand)
        first_real_integral = self._integrator.integrate(integrand)
        return first_real_integral

    def evaluate_second_real_integral(
        self, gap_energy: float, temperature: float, omega: float
    ) -> float:
        if h_bar * omega <= 2 * gap_energy:
            return 0

        integrand = MattisBardeenRealSecondIntegrand(gap_energy, temperature, omega)
        integrand = remove_lower_singularity(integrand)
        integrand = remove_upper_singularity(integrand)
        second_real_integral = self._integrator.integrate(integrand)
        return second_real_integral

    def evaluate_imaginary_integral(
        self, gap_energy: float, temperature: float, omega: float
    ) -> float:

        integrand = MattisBardeenImaginaryIntegrand(gap_energy, temperature, omega)
        integrand = remove_lower_singularity(integrand)
        integrand = remove_upper_singularity(integrand)
        imaginary_integral = self._integrator.integrate(integrand)
        return imaginary_integral

    def evaluate(self, temperature: float, frequency: float) -> complex:
        omega = 2 * pi * frequency
        gap_energy = self._gap_energy.evaluate(temperature)

        sigma_r1 = self.evaluate_first_real_integral(gap_energy, temperature, omega)
        sigma_r2 = self.evaluate_second_real_integral(gap_energy, temperature, omega)
        sigma_i = self.evaluate_imaginary_integral(gap_energy, temperature, omega)

        scale = self._conductivity_0 / (h_bar * omega)
        unscaled = 2 * sigma_r1 + sigma_r2 + sigma_i * 1j

        return scale * unscaled


__all__ = ["MattisBardeenComplexConductivity"]
