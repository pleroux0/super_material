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

    Solves the self consistent gap energy equation :cite:`BCSTheory`

    .. math::
        \\frac
            {1}
            {N(0) V}
        =
        \\int\\limits_{0}^{\\hbar \\omega_{D}}
            \\frac
                {\\tanh{
                    \\left(
                    \\frac{\\sqrt{\\Delta^{2} + x^{2}}}{2 T k_{B}}
                    \\right)}}
                {\\sqrt{\\Delta^{2} + x^{2}}}
            dx

    The Debye frequency is given by :math:`\\omega_D`. For convenience we
    define :math:`\\kappa = \\frac{\\hbar \\omega_{D}}{2 T_{c} k_{B}}` and
    :math:`\\eta = \\frac{1}{N(0) V}`.

    The original BCS assumes the weak coupling limit which corresponds to
    :math:`\\kappa \\gg 1`. We do not make the assumption when calculting the
    gap energy.
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
        temperature_scale = 1 / (2 * k_B * temperature)

        def gap_energy_integrand(z, dirac):
            expr = sqrt(dirac * dirac + z * z)
            return tanh(expr * temperature_scale) / expr

        def equation(dirac):
            result = integrate.quad(gap_energy_integrand, 0, to, args=(dirac))
            return self.eta() - result[0]

        result = optimize.bisect(equation, tol, self.critical_temperature() - tol)

        return abs(result)

    def critical_temperature(self) -> float:
        return self.gap_energy_0() * sinh(self.eta()) / (2 * self.kappa() * k_B)


__all__ = ["BCSGapEnergy"]
