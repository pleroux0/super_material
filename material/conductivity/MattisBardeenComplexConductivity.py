from math import sqrt, exp

from ..integrate import IntegrandBoundary, IntegrandInterface, IntegrandInterval

from ..constants import h_bar, k_B


def f(E, T):
    exponent = E / (k_B * T)

    if exponent > 500:
        return 0

    denumerator = exp(exponent) + 1
    return 1 / denumerator


class MattisBardeenImaginaryIntegrand(IntegrandInterface):
    _gap_energy_0: float
    _temperature: float
    _omega: float

    def __init__(self, gap_energy_0: float, temperature: float, omega: float):
        self._gap_energy_0 = gap_energy_0
        self._temperature = temperature
        self._omega = omega

    def interval(self) -> IntegrandInterval:
        lower = self._gap_energy_0 - h_bar * self._omega

        if lower > -self._gap_energy_0:
            start = IntegrandBoundary(lower, False)
        else:
            start = IntegrandBoundary(-self._gap_energy_0, True)

        end = IntegrandBoundary(self._gap_energy_0, False)

        return IntegrandInterval(start, end)

    def evaluate(self, E: float) -> float:
        a = 1 - 2 * f(E + h_bar * self._omega, self._temperature)
        b = E ** 2 + self._gap_energy_0 ** 2 + h_bar * self._omega * E
        c = sqrt(self._gap_energy_0 ** 2 - E ** 2)
        d = sqrt((E + h_bar * self._omega) ** 2 - self._gap_energy_0 ** 2)

        return a * b / (c * d)
