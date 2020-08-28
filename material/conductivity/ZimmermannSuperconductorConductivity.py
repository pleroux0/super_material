from math import sqrt, exp, pi, inf, sin, tanh
from typing import List
from dataclasses import dataclass

import scipy.integrate
import numpy as np

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


@dataclass
class AbsOrRelTolerance:
    atol: float = 1e-6
    rtol: float = 1e-6

    def __call__(self, previous, current):
        return (
            abs(previous - current) < self.rtol * previous
            or abs(previous - current) < self.atol
        )


class ChebychevGaussQuadrature:
    _x: List[np.ndarray]
    _w: List[np.ndarray]
    _num_levels: int

    __slots__ = ("_x", "_w", "_num_levels")

    def __init__(self):
        self._x = []
        self._w = []
        self._num_levels = 0

    def _next_level(self):
        x, w = scipy.special.roots_chebyt(2 ** self._num_levels)
        self._x.append(x)
        self._w.append(w)
        self._num_levels += 1

    def _get_level(self, level: int):
        assert level >= 0

        while level >= self._num_levels:
            self._next_level()

        return self._x[level], self._w[level]

    def _integrate_level(self, integrand, level):
        x, w = self._get_level(level)
        return integrand(x).dot(w)

    def integrate(
        self,
        f,
        a=-1,
        b=1,
        min_levels=10,
        max_levels=20,
        converge=AbsOrRelTolerance(1e-6),
    ):
        m = (b - a) * 0.5
        c = (b + a) * 0.5

        def integrand(x):
            return m * f(m * x + c)

        prev_result = self._integrate_level(integrand, 0)
        level = 1

        while level <= max_levels:
            result = self._integrate_level(integrand, level)

            if level >= min_levels and converge(prev_result, result):
                return result

            prev_result = result
            level += 1

        raise Exception(
            "Result did not converge: prev_result={} result={}".format(
                prev_result, result
            )
        )


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


def _p1(delta, e, w):
    return np.sqrt((e + w * h_bar) ** 2 - delta ** 2)


def _p2(delta, e):
    return np.sqrt(e ** 2 - delta ** 2)


def _p3(delta, e, w):
    return np.sqrt((e - w * h_bar) ** 2 - delta ** 2)


def _p4(delta, e, w):
    return np.sqrt(delta ** 2 - (e - w * h_bar) ** 2)


def _th1(e, T):
    return np.tanh(e / (2 * k_B * T))


def _th2(e, w, T):
    return np.tanh((w * h_bar + e) / (2 * k_B * T))


def _integrand_i1(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p4 = 1j * _p4(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau

    c42 = (delta ** 2 + e * (e - w * h_bar)) / (p2 * p4)

    return th * ((1 - c42) / (p4 + p2 + itau) - (1 + c42) / (p4 - p2 + itau))


def _integrand_i1_konrad(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p4 = 1j * _p4(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau

    return th * ((1) / (p4 + p2 + itau) - (1) / (p4 - p2 + itau))


def _integrand_i1_cheby_part_1(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p4 = 1j * _p4(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau

    c42 = (
        (delta ** 2 + e * (e - w * h_bar))
        / (1j * np.sqrt((e + delta) * (delta + e - h_bar * w)))
        * 2
        / (w * h_bar)
    )

    return th * ((-c42) / (p4 + p2 + itau) - (+c42) / (p4 - p2 + itau))


def _integrand_i1_cheby_part_2(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p4 = 1j * _p4(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau

    c42 = (delta ** 2 + e * (e - w * h_bar)) / (1j * p2 * delta)

    return th * ((-c42) / (p4 + p2 + itau) - (+c42) / (p4 - p2 + itau))


def _integrand_i3(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p3 = _p3(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau
    c32 = (delta ** 2 + e * (e - w * h_bar)) / (p2 * p3)

    return th * ((1 - c32) / (p3 + p2 + itau) - (1 + c32) / (p3 - p2 + itau))


def _integrand_i3_konrad(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p3 = _p3(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau

    return th * ((1) / (p3 + p2 + itau) - (1) / (p3 - p2 + itau))


def _integrand_i3_cheby(e, delta, w, tau, T):
    p2 = _p2(delta, e)
    p3 = _p3(delta, e, w)
    th = _th1(e, T)

    itau = 1j / tau
    c32 = (
        (delta ** 2 + e * (e - w * h_bar))
        / (np.sqrt(-(delta + e) * (e - delta - h_bar * w)))
        * 2
        / (w * h_bar - 2 * delta)
    )

    return th * ((-c32) / (p3 + p2 + itau) - (c32) / (p3 - p2 + itau))


def _integrand_i2(e, delta, w, tau, T):
    p1 = _p1(delta, e, w)
    p2 = _p2(delta, e)

    th1 = _th1(e, T)
    th2 = _th2(e, w, T)

    itau = 1j / tau
    c12 = (delta ** 2 + e * (e + w * h_bar)) / (p1 * p2)

    out = th2 * ((1 + c12) / (p1 - p2 + itau) - (1 - c12) / (-p1 - p2 + itau))
    out += th1 * ((1 - c12) / (p1 + p2 + itau) - (1 + c12) / (p1 - p2 + itau))

    return out


class ZimmermansConductivity:
    _gap_energy: GapEnergyInterface
    _tau: float
    _conductivity_0: float
    _chebyshev_gauss_quadrature: ChebychevGaussQuadrature

    __slots__ = (
        "_gap_energy",
        "_tau",
        "_conductivity_0",
        "_chebyshev_gauss_quadrature",
    )

    def __init__(
        self, gap_energy: GapEnergyInterface, tau: float, conductivity_0: float
    ):
        self._gap_energy = gap_energy
        self._tau = tau / h_bar
        self._conductivity_0 = conductivity_0
        self._chebyshev_gauss_quadrature = ChebychevGaussQuadrature()

    def evaluate(self, frequency: float = 10e9, temperature: float = 4.2) -> complex:
        delta = self._gap_energy.evaluate(temperature)
        w = frequency * 2 * np.pi
        args = (w, self._tau, delta, temperature)

        tol = 1.49e-08
        # tol = 1e-12

        s = self._calc_i2(*args, tol)

        if h_bar * w <= 2 * delta:
            s += self._calc_i1_part_1(*args, tol)
        else:
            s += self._calc_i1_part_2(*args, tol)
            s += self._calc_i3(*args, tol)

        return s * 1j / (2 * h_bar * w * self._tau) * self._conductivity_0

    def _calc_i1_part_1(self, w, tau, gap_energy, T, tol=1e-6):
        delta = gap_energy

        a, b = delta, delta + h_bar * w
        args = (delta, w, tau, T)

        r_konrad, _ = scipy.integrate.quadrature(
            _integrand_i1_konrad,
            a,
            b,
            args=args,
            tol=tol,
            rtol=tol,
            miniter=10,
            maxiter=200,
        )

        # r_konrad, _ = scipy.integrate.quad(
        #     lambda x: _integrand_i1_konrad(x, *args),
        #     a,
        #     b,
        #     epsabs=tol,
        #     epsrel=tol,
        #     points=(a, b),
        # )

        def integrand_cheby(e):
            return _integrand_i1_cheby_part_1(e, *args)

        r_cheby = self._chebyshev_gauss_quadrature.integrate(
            integrand_cheby, a, b, min_levels=10, max_levels=20
        )

        return r_cheby + r_konrad

    def _calc_i1_part_2(self, w, tau, gap_energy, T, tol=1e-6):
        delta = gap_energy

        a, b = h_bar * w - delta, delta + h_bar * w
        args = (delta, w, tau, T)

        r_konrad, _ = scipy.integrate.quadrature(
            _integrand_i1_konrad, a, b, args=args, tol=tol, miniter=20, maxiter=200
        )

        def integrand_cheby(e):
            return _integrand_i1_cheby_part_2(e, *args)

        r_cheby = self._chebyshev_gauss_quadrature.integrate(
            integrand_cheby, a, b, min_levels=10, max_levels=20
        )

        return r_cheby + r_konrad

    def _calc_i3(self, w, tau, gap_energy, T, tol=1e-6):
        delta = gap_energy

        a, b = delta, -delta + h_bar * w
        args = (delta, w, tau, T)

        r_konrad, _ = scipy.integrate.quadrature(
            _integrand_i3_konrad,
            a,
            b,
            args=args,
            tol=tol,
            rtol=tol,
            miniter=20,
            maxiter=200,
        )

        def integrand_cheby(e):
            return _integrand_i3_cheby(e, *args)

        r_cheby = self._chebyshev_gauss_quadrature.integrate(
            integrand_cheby, a, b, min_levels=10, max_levels=20
        )

        return r_cheby + r_konrad

    @staticmethod
    def _calc_i2(w, tau, gap_energy, T, tol=1e-6):
        delta = gap_energy

        def integrand(x):
            m = delta + (x / (1 - x)) ** 2
            suffix = 2 * x / ((1 - x) ** 3)
            return _integrand_i2(m, delta, w, tau, T) * suffix

        out = scipy.integrate.quadrature(
            integrand, 0, 1, tol=tol, miniter=20, maxiter=200, rtol=tol
        )

        return out[0]


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
            self._gap_energy, self._scattering_time, self._temperature, self._omega,
        )
        return equations.I2(E)


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
        upper = IntegrandBoundary(1, False)
        interval = IntegrandInterval(lower, upper)
        return interval

    def evaluate(self, x: float) -> float:
        E = self._gap_energy + (x / (1 - x)) ** 2
        scale = 2 * x / ((1 - x) ** 3)
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
    _gap_energy: GapEnergyInterface
    _conductivity_0: float
    _scattering_time: float
    _integrator: ScipyQuadratureIntegrator

    _legacy: ZimmermansConductivity

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
            absolute_tolerance=1e-8, relative_tolerance=1e-6, maximum_order=200
        )
        self._legacy = ZimmermansConductivity(
            gap_energy, scattering_time, conductivity_0
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
        args = (omega, self._scattering_time / h_bar, gap_energy, temperature)
        tol = 1.49e-08
        s = self._legacy._calc_i2(*args, tol)

        return s

        integrand = ZimmermannSecondIntegralTransformed(
            gap_energy, self._scattering_time, temperature, omega
        )
        integrand = remove_lower_singularity(integrand)
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
        s = 0

        args = (omega, self._scattering_time / h_bar, gap_energy, temperature)
        tol = 1.49e-08

        if h_bar * omega <= 2 * gap_energy:
            s += self._legacy._calc_i1_part_1(*args, tol)
        else:
            s += self._legacy._calc_i1_part_2(*args, tol)
            s += self._legacy._calc_i3(*args, tol)

        return s

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

        scale = self._conductivity_0 * 1j / (2 * omega * self._scattering_time)

        J = self.evaluate_j(gap_energy, temperature, omega)
        second_integral = self.evaluate_second_integral(gap_energy, temperature, omega)

        out = scale * (J + second_integral)

        return out
