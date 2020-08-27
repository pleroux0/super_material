from math import sqrt, exp, pi, inf, sin, tanh

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
)

from ..constants import h_bar, k_B


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

        tmp = self._gap_energy ** 2 + E * (E - self._omega * h_bar) / (p2 * p4)

        out = th1 * (
            (1 - tmp) / (p4 + p2 + 1j / self._scattering_time)
            - (1 + tmp) / (p4 - p2 + 1j / self._scattering_time)
        )

        return out

    def I2(self, E: float) -> complex:
        p1 = self.p1(E)
        p2 = self.p2(E)
        th1 = self.th1(E)
        th2 = self.th2(E)

        tmp = self._gap_energy ** 2 + E * (E + self._omega * h_bar) / (p1 * p2)

        part1 = th2 * (
            (1 + tmp) / (p1 - p2 + 1j / self._scattering_time)
            - (1 - tmp) / (-p1 - p2 + 1j / self._scattering_time)
        )

        part2 = th1 * (
            (1 - tmp) / (p1 + p2 + 1j / self._scattering_time)
            - (1 + tmp) / (p1 - p2 + 1j / self._scattering_time)
        )

        return part1 + part2

    def I3(self, E: float) -> complex:
        p2 = self.p2(E)
        p3 = self.p3(E)
        th1 = self.th1(E)

        tmp = self._gap_energy ** 2 + E * (E - self._omega * h_bar) / (p3 * p2)

        out = th1 * (
            (1 - tmp) / (p3 + p2 + 1j / self._scattering_time)
            - (1 + tmp) / (p3 - p2 + 1j / self._scattering_time)
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
            self._gap_energy, self._temperature, self._omega, self._scattering_time
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
            self._gap_energy, self._temperature, self._omega, self._scattering_time
        )
        return equations.I1(E)


class ZimmermannSecondIntegral(IntegrandInterface):
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
        upper = IntegrandBoundary(inf, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, E: float) -> float:
        equations = ZimmermannSuperconductorEquations(
            self._gap_energy, self._temperature, self._omega, self._scattering_time
        )
        return equations.I2(E)


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
            self._gap_energy, self._temperature, self._omega, self._scattering_time
        )
        return equations.I3(E)


class ZimmermannSuperconductorConductivity(SuperconductorConductivityInterface):
    _gap_energy: GapEnergyInterface
    _conductivity_0: float
    _scattering_time: float
    _integrator: QuadpackIntegrator

    def __init__(
        self,
        gap_energy: GapEnergyInterface,
        conductivity_0: float,
        scattering_time: float,
    ):
        self._gap_energy = gap_energy
        self._conductivity_0 = conductivity_0
        self._scattering_time = scattering_time
        self._integrator = QuadpackIntegrator(
            absolute_tolerance=1e-8, relative_tolerance=1e-8, limit=50
        )

    def evaluate_first_integral_superconductor_part(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannFirstIntegralSuperconductorPart(
            gap_energy, self._scattering_time, temperature, omega
        )
        first_integral = self._integrator.integrate(integrand)
        return first_integral

    def evaluate_first_integral_normal_part(
        self, gap_energy: float, temperature: float, omega: float
    ):
        integrand = ZimmermannFirstIntegralNormalPart(
            gap_energy, self._scattering_time, temperature, omega
        )
        first_integral = self._integrator.integrate(integrand)
        return first_integral

    def evaluate_second_integral(self, gap_energy: float, temperature: float, omega: float):
        integrand = ZimmermannSecondIntegral(
            gap_energy, self._scattering_time, temperature, omega
        )
        second_integral = self._integrator.integrate(integrand)
        return second_integral

    def evaluate_third_integral(self, gap_energy: float, temperature: float, omega: float):
        integrand = ZimmermannThirdIntegral(
            gap_energy, self._scattering_time, temperature, omega
        )
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

        scale = self._conductivity_0 * 1j / (2 * h_bar * omega * self._scattering_time)

        J = self.evaluate_j(gap_energy, temperature, omega)
        second_integral = self.evaluate_second_integral(gap_energy, temperature, omega)

        out = scale * (J + second_integral)

        return out
