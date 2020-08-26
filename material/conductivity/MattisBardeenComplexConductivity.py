from math import sqrt, exp, pi

from .SuperconductorConductivityInterface import SuperconductorConductivityInterface

from ..gap_energy.GapEnergyInterface import GapEnergyInterface
from ..integrate import (
    IntegrandBoundary,
    IntegrandInterface,
    IntegrandInterval,
    IntegrandIntervalTransformInterface,
    QuadpackIntegrator,
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


class MattisBardeenRealFirstIntegrand(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float

    def __init__(self, gap_energy: float, temperature: float, omega: float):
        assert h_bar * omega > 2 * gap_energy

        self._gap_energy = gap_energy
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = self._gap_energy - h_bar * self._omega
        upper = -self._gap_energy
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
        lower = self._gap_energy - h_bar * self._omega
        upper = -self._gap_energy
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

    def __init__(self, gap_energy: GapEnergyInterface, conductivity_0: float):
        self._gap_energy = gap_energy
        self._conductivity_0 = conductivity_0

    @staticmethod
    def evaluate_first_real_integral(
        gap_energy: float, temperature: float, omega: float
    ) -> float:
        integrator = QuadpackIntegrator()
        integrand = MattisBardeenRealFirstIntegrand(gap_energy, temperature, omega)
        second_real_integral = integrator.integrate(integrand)
        return second_real_integral

    @staticmethod
    def evaluate_second_real_integral(
        gap_energy: float, temperature: float, omega: float
    ) -> float:
        integrator = QuadpackIntegrator()
        integrand = MattisBardeenRealSecondIntegrand(gap_energy, temperature, omega)
        first_real_integral = integrator.integrate(integrand)
        return first_real_integral

    @staticmethod
    def evaluate_imaginary_integral(
        gap_energy: float, temperature: float, omega: float
    ) -> float:
        integrator = QuadpackIntegrator()
        integrand = MattisBardeenImaginaryIntegrand(gap_energy, temperature, omega)
        imaginary_integral = integrator.integrate(integrand)
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
