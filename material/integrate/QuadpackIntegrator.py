from scipy.integrate import quad

from .IntegratorInterface import IntegratorInterface
from .IntegrandInterface import IntegrandInterface


class QuadpackIntegrator(IntegratorInterface):
    def integrate(self, integrand: IntegrandInterface) -> float:
        interval = integrand.interval()
        start = interval.start().value()
        end = interval.end().value()
        f = integrand.evaluate

        output, _ = quad(f, start, end)

        return output


__all__ = ["QuadpackIntegrator"]
