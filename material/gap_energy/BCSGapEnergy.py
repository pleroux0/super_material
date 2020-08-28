from math import tanh, sinh, sqrt
from functools import lru_cache

import scipy.integrate as integrate
import scipy.optimize as optimize

from ..constants import k_B
from ..integrate.IntegrandBoundary import IntegrandBoundary
from ..integrate.IntegrandInterface import IntegrandInterface
from ..integrate.IntegrandInterval import IntegrandInterval
from ..integrate.QuadpackIntegrator import QuadpackIntegrator

from .GapEnergyInterface import GapEnergyInterface


class BCSEtaIntegrand(IntegrandInterface):
    _kappa: float

    def __init__(self, kappa):
        self._kappa = kappa

    def kappa(self) -> float:
        return self._kappa

    def evaluate(self, x: float) -> float:
        return tanh(x) / x

    def interval(self) -> IntegrandInterval:
        start = IntegrandBoundary(0, False)
        end = IntegrandBoundary(self.kappa(), True)
        interval = IntegrandInterval(start, end)
        return interval


class BCSGapEnergy(GapEnergyInterface):
    """ Gap energy as calcuated from BCS theory

    Solves the self consistent gap energy equation

    """

    _gap_energy_0: float  # In Electron Volt
    _kappa: float

    # No dynamic variables
    __slots__ = ()

    def __init__(self, gap_energy_0: float, kappa: float):
        assert gap_energy_0 > 0
        assert kappa > 0

        self._gap_energy_0 = gap_energy_0
        self._kappa = kappa

    def gap_energy_0(self):
        return self._gap_energy_0

    def kappa(self):
        return self._kappa

    @lru_cache
    def eta(self):
        """ Calculate and return eta """

        integrand = BCSEtaIntegrand(self.kappa())
        integrator = QuadpackIntegrator()
        eta = integrator.integrate(integrand)
        return eta

    def evaluate(self, temperature: float) -> float:
        assert temperature >= 0
        tol = 1e-12

        if temperature < tol:
            return self.gap_energy_0()

        T_c = self.critical_temperature()
        if temperature > T_c - tol:
            return 0

        to = self.gap_energy_0() * sinh(self.eta())

        def gap_energy_integrand(z, dirac):
            dirac2 = dirac * dirac
            z2 = z * z
            expr = sqrt(dirac2 + z2)
            return tanh(expr / (2 * k_B * temperature)) / expr

        def equation(dirac):
            result = integrate.quad(gap_energy_integrand, 0, to, args=(dirac))
            return self.eta() - result[0]

        result = optimize.bisect(equation, tol, self.critical_temperature() - tol)

        return abs(result)

    def critical_temperature(self) -> float:
        return self.gap_energy_0() * sinh(self.eta()) / (2 * self.kappa() * k_B)


__all__ = ["BCSGapEnergy"]
