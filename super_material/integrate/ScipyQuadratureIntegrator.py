from scipy.integrate import quadrature
from numpy import vectorize

from .IntegratorInterface import IntegratorInterface
from .IntegrandInterface import IntegrandInterface


class ScipyQuadratureIntegrator(IntegratorInterface):
    _absolute_tolerance: float
    _relative_tolerance: float
    _maximum_order: int

    def __init__(
        self,
        absolute_tolerance: float = 1.49e-8,
        relative_tolerance: float = 1.49e-8,
        maximum_order: int = 50,
        minimum_order: int = 1,
    ):
        self._absolute_tolerance = absolute_tolerance
        self._relative_tolerance = relative_tolerance
        self._maximum_order = maximum_order
        self._minimum_order = minimum_order

    def integrate(self, integrand: IntegrandInterface) -> float:
        interval = integrand.interval()
        start = interval.start().value()
        end = interval.end().value()
        f = integrand.evaluate

        output, _ = quadrature(
            vectorize(f),
            start,
            end,
            tol=self._absolute_tolerance,
            rtol=self._relative_tolerance,
            maxiter=self._maximum_order,
            miniter=self._minimum_order,
        )

        return output


__all__ = ["ScipyQuadratureIntegrator"]
