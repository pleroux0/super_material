from math import sqrt, exp, pi, inf, sin

from .SuperconductorConductivityInterface import SuperconductorConductivityInterface

from ..gap_energy.GapEnergyInterface import GapEnergyInterface
from ..integrate import (
    IntegrandBoundary,
    IntegrandInterface,
    IntegrandInterval,
    QuadpackIntegrator,
    TransformedIntegrand,
    ChebyshevLowerSingularityTransform,
    ChebyshevUpperSingularityTransform,
    ChebyshevSingularityTransform,
)

from ..constants import h_bar, k_B


def fermi_dirac_function(E, T):
    exponent = E / (k_B * T)

    if exponent > 500:
        return 0

    denumerator = exp(exponent) + 1
    return 1 / denumerator


def remove_chebyshev_singularity(integrand: IntegrandInterface):
    lower_singularity = integrand.interval().start().value()
    upper_singularity = integrand.interval().end().value()
    transform = ChebyshevSingularityTransform(lower_singularity, upper_singularity)
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
        lower = IntegrandBoundary(0, False)
        upper = IntegrandBoundary(1 / (self._gap_energy ** 2), True)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, x: float) -> float:
        E = self._gap_energy / sqrt(1 - (self._gap_energy ** 4) * (x ** 2))

        a = fermi_dirac_function(E, self._temperature) - fermi_dirac_function(
            E + h_bar * self._omega, self._temperature
        )
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy ** 2)

        return E ** 2 * a * b / d


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
        lower = IntegrandBoundary(0, True)
        upper = IntegrandBoundary(pi, True)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, x: float) -> float:
        assert 0 <= x <= pi

        lower = self._gap_energy - h_bar * self._omega
        upper = -self._gap_energy
        E = lower + (upper - lower) * (sin(x / 2) ** 2)

        a = 1 - 2 * fermi_dirac_function(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        c = sqrt(-E + self._gap_energy)
        d = sqrt((E + h_bar * self._omega) + self._gap_energy)

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
        a = 1 - 2 * fermi_dirac_function(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy ** 2 + h_bar * self._omega * E
        c = sqrt(self._gap_energy ** 2 - E ** 2)
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy ** 2)

        return a * b / (c * d)


class MattisBardeenSuperconductorConductivity(SuperconductorConductivityInterface):
    """ Superconductor conductivity as calculated by Mattis and Bardeen

    Numerically evaluates the integral expression of Mattis and Bardeen
    :cite:`MattisBardeenSuperconductorConductivity`
    """

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
        # return 0
        integrand = MattisBardeenRealFirstIntegrand(gap_energy, temperature, omega)
        first_real_integral = self._integrator.integrate(integrand)
        return first_real_integral

    def evaluate_second_real_integral(
        self, gap_energy: float, temperature: float, omega: float
    ) -> float:
        if h_bar * omega <= 2 * gap_energy:
            return 0

        integrand = MattisBardeenRealSecondIntegrand(gap_energy, temperature, omega)
        second_real_integral = self._integrator.integrate(integrand)
        return second_real_integral

    def evaluate_imaginary_integral(
        self, gap_energy: float, temperature: float, omega: float
    ) -> float:

        integrand = MattisBardeenImaginaryIntegrand(gap_energy, temperature, omega)
        integrand = remove_chebyshev_singularity(integrand)
        imaginary_integral = self._integrator.integrate(integrand)
        return imaginary_integral

    def evaluate(self, temperature: float, frequency: float) -> complex:
        omega = 2 * pi * frequency
        gap_energy = self._gap_energy.evaluate(temperature)

        sigma_r1 = self.evaluate_first_real_integral(gap_energy, temperature, omega)
        sigma_r2 = self.evaluate_second_real_integral(gap_energy, temperature, omega)
        sigma_i = self.evaluate_imaginary_integral(gap_energy, temperature, omega)

        scale = self._conductivity_0 / (h_bar * omega)
        unscaled = 2 * sigma_r1 - sigma_r2 + sigma_i * 1j

        return scale * unscaled


__all__ = ["MattisBardeenComplexConductivity"]
