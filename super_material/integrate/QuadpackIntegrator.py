from scipy.integrate import quad

from .IntegratorInterface import IntegratorInterface
from .IntegrandInterface import IntegrandInterface


class QuadpackIntegrator(IntegratorInterface):
    _absolute_tolerance: float
    _relative_tolerance: float
    _limit: int

    def __init__(
        self,
        absolute_tolerance: float = 1.49e-8,
        relative_tolerance: float = 1.49e-8,
        limit: int = 50,
    ):
        self._absolute_tolerance = absolute_tolerance
        self._relative_tolerance = relative_tolerance
        self._limit = limit

    def integrate(self, integrand: IntegrandInterface) -> float:
        interval = integrand.interval()
        start = interval.start().value()
        end = interval.end().value()
        f = integrand.evaluate

        output, _ = quad(
            f,
            start,
            end,
            epsabs=self._absolute_tolerance,
            epsrel=self._relative_tolerance,
            limit=self._limit,
        )

        return output


__all__ = ["QuadpackIntegrator"]
