from math import sqrt, pi, tanh, cos

from .SuperconductorConductivityInterface import SuperconductorConductivityInterface

from ..gap_energy.GapEnergyInterface import GapEnergyInterface
from ..integrate import (
    IntegrandBoundary,
    IntegrandInterface,
    IntegrandInterval,
    ScipyQuadratureIntegrator,
    TransformedIntegrand,
    IntegrandIntervalTransformInterface,
    ChebyshevLowerSingularityTransform,
    ChebyshevUpperSingularityTransform,
)

from ..constants import h_bar, k_B


def remove_lower_singularity(integrand: IntegrandInterface):
    lower_singularity = integrand.interval().start().value()
    transform = ChebyshevLowerSingularityTransform(lower_singularity)
    transformed_integrand = TransformedIntegrand(integrand, transform)
    return transformed_integrand


def remove_upper_singularity(integrand: IntegrandInterface):
    upper_singularity = integrand.interval().end().value()
    transform = ChebyshevUpperSingularityTransform(upper_singularity)
    transformed_integrand = TransformedIntegrand(integrand, transform)
    return transformed_integrand


class ZimmermannSuperconductorEquations:
    _gap_energy: float
    _temperature: float
    _omega: float
    _scattering_time: float

    def __init__(
        self,
        gap_energy: float,
        scattering_time: float,
        temperature: float,
        omega: float,
    ):
        self._gap_energy = gap_energy
        self._scattering_time = scattering_time
        self._temperature = temperature
        self._omega = omega

    def p1(self, E: float) -> float:
        return sqrt((E + self._omega * h_bar) ** 2 - self._gap_energy ** 2)

    def p2(self, E: float) -> float:
        return sqrt(E ** 2 - self._gap_energy ** 2)

    def p3(self, E: float) -> float:
        return sqrt((E - self._omega * h_bar) ** 2 - self._gap_energy ** 2)

    def p4(self, E: float) -> complex:
        return 1j * sqrt(self._gap_energy ** 2 - (E - self._omega * h_bar) ** 2)

    def th1(self, E: float) -> float:
        return tanh(E / (2 * k_B * self._temperature))

    def th2(self, E: float) -> float:
        return tanh((E + self._omega * h_bar) / (2 * k_B * self._temperature))

    def I1(self, E: float) -> complex:
        p2 = self.p2(E)
        p4 = self.p4(E)
        th1 = self.th1(E)

        tmp = (self._gap_energy ** 2 + E * (E - self._omega * h_bar)) / (p2 * p4)

        out = th1 * (
            (1 - tmp) / ((p4 + p2) * self._scattering_time + h_bar * 1j)
            - (1 + tmp) / ((p4 - p2) * self._scattering_time + h_bar * 1j)
        )

        return out

    def I2(self, E: float) -> complex:
        p1 = self.p1(E)
        p2 = self.p2(E)
        th1 = self.th1(E)
        th2 = self.th2(E)

        tmp = (self._gap_energy ** 2 + E * (E + self._omega * h_bar)) / (p1 * p2)

        part1 = th2 * (
            (1 + tmp) / ((p1 - p2) * self._scattering_time + h_bar * 1j)
            - (1 - tmp) / ((-p1 - p2) * self._scattering_time + h_bar * 1j)
        )

        part2 = th1 * (
            (1 - tmp) / ((p1 + p2) * self._scattering_time + h_bar * 1j)
            - (1 + tmp) / ((p1 - p2) * self._scattering_time + h_bar * 1j)
        )

        return part1 + part2

    def I3(self, E: float) -> complex:
        p2 = self.p2(E)
        p3 = self.p3(E)
        th1 = self.th1(E)

        tmp = (self._gap_energy ** 2 + E * (E - self._omega * h_bar)) / (p3 * p2)

        out = th1 * (
            (1 - tmp) / ((p3 + p2) * self._scattering_time + h_bar * 1j)
            - (1 + tmp) / ((p3 - p2) * self._scattering_time + h_bar * 1j)
        )

        return out


class ZimmermannFirstIntegralSuperconductorPart(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float
    _scattering_time: float

    __slots__ = ()

    def __init__(
        self,
        gap_energy: float,
        scattering_time: float,
        temperature: float,
        omega: float,
    ):
        self._gap_energy = gap_energy
        self._scattering_time = scattering_time
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(self._gap_energy, False)
        upper = IntegrandBoundary(h_bar * self._omega + self._gap_energy, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        equations = ZimmermannSuperconductorEquations(
            self._gap_energy, self._scattering_time, self._temperature, self._omega,
        )
        return equations.I1(E)


class ZimmermannFirstIntegralNormalPart(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float
    _scattering_time: float

    __slots__ = ()

    def __init__(
        self,
        gap_energy: float,
        scattering_time: float,
        temperature: float,
        omega: float,
    ):
        self._gap_energy = gap_energy
        self._scattering_time = scattering_time
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(h_bar * self._omega - self._gap_energy, False)
        upper = IntegrandBoundary(h_bar * self._omega + self._gap_energy, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        equations = ZimmermannSuperconductorEquations(
            self._gap_energy, self._scattering_time, self._temperature, self._omega,
        )
        return equations.I1(E)


class ZimmermannSecondIntegralTransformed(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float
    _scattering_time: float

    __slots__ = ()

    def __init__(
        self,
        gap_energy: float,
        scattering_time: float,
        temperature: float,
        omega: float,
    ):
        self._gap_energy = gap_energy
        self._scattering_time = scattering_time
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(0, False)
        upper = IntegrandBoundary(pi / (2 * self._gap_energy), False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, x: float) -> float:
        E = self._gap_energy / cos(self._gap_energy * x)
        scale = E * sqrt(E ** 2 - self._gap_energy ** 2)
        equations = ZimmermannSuperconductorEquations(
            self._gap_energy, self._scattering_time, self._temperature, self._omega,
        )
        I2 = equations.I2(E)
        return I2 * scale


class ZimmermannThirdIntegral(IntegrandInterface):
    _gap_energy: float
    _temperature: float
    _omega: float
    _scattering_time: float

    __slots__ = ()

    def __init__(
        self,
        gap_energy: float,
        scattering_time: float,
        temperature: float,
        omega: float,
    ):
        self._gap_energy = gap_energy
        self._scattering_time = scattering_time
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = IntegrandBoundary(self._gap_energy, False)
        upper = IntegrandBoundary(h_bar * self._omega - self._gap_energy, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        equations = ZimmermannSuperconductorEquations(
            self._gap_energy, self._scattering_time, self._temperature, self._omega,
        )
        return equations.I3(E)


class ZimmermannSuperconductorConductivity(SuperconductorConductivityInterface):
    """ Superconductor conductivity as calculated by Zimmermann

    Numerically evaluates the integral expression of Zimmermann
    :cite:`ZimmermannSuperconductorConductivity`
    """

    _gap_energy: GapEnergyInterface
    _conductivity_0: float
    _scattering_time: float
    _integrator: ScipyQuadratureIntegrator

    def __init__(
        self,
        gap_energy: GapEnergyInterface,
        conductivity_0: float,
        scattering_time: float,
    ):
        self._gap_energy = gap_energy
        self._conductivity_0 = conductivity_0
        self._scattering_time = scattering_time
        self._integrator = ScipyQuadratureIntegrator(
            absolute_tolerance=1e-6,
            relative_tolerance=1e-6,
            maximum_order=200,
            minimum_order=1,
        )

    def evaluate_first_integral_superconductor_part(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannFirstIntegralSuperconductorPart(
            gap_energy, self._scattering_time, temperature, omega
        )
        integrand = remove_lower_singularity(integrand)
        integrand = remove_upper_singularity(integrand)
        first_integral = self._integrator.integrate(integrand)
        return first_integral

    def evaluate_first_integral_normal_part(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannFirstIntegralNormalPart(
            gap_energy, self._scattering_time, temperature, omega
        )
        integrand = remove_lower_singularity(integrand)
        integrand = remove_upper_singularity(integrand)
        first_integral = self._integrator.integrate(integrand)
        return first_integral

    def evaluate_second_integral(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannSecondIntegralTransformed(
            gap_energy, self._scattering_time, temperature, omega
        )
        second_integral = self._integrator.integrate(integrand)
        return second_integral

    def evaluate_third_integral(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannThirdIntegral(
            gap_energy, self._scattering_time, temperature, omega
        )
        integrand = remove_lower_singularity(integrand)
        integrand = remove_upper_singularity(integrand)
        third_integral = self._integrator.integrate(integrand)
        return third_integral

    def evaluate_j(self, gap_energy: float, temperature: float, omega: float):
        if h_bar * omega <= 2 * gap_energy:
            first_integral = self.evaluate_first_integral_superconductor_part(
                gap_energy, temperature, omega
            )
            return first_integral

        third_integral = self.evaluate_third_integral(gap_energy, temperature, omega)
        first_integral = self.evaluate_first_integral_normal_part(
            gap_energy, temperature, omega
        )

        J = third_integral + first_integral

        return J

    def evaluate(self, temperature: float, frequency: float) -> complex:
        omega = 2 * pi * frequency
        gap_energy = self._gap_energy.evaluate(temperature)

        scale = self._conductivity_0 * 1j / (2 * omega)

        J = self.evaluate_j(gap_energy, temperature, omega)
        second_integral = self.evaluate_second_integral(gap_energy, temperature, omega)

        out = scale * (J + second_integral)

        return out
