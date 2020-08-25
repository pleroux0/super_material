from typing import Callable

from .IntegrandInterface import IntegrandInterface
from .IntegrandInterval import IntegrandInterval


class BaseIntegrand(IntegrandInterface):
    _function: Callable[[float], float]
    _interval: IntegrandInterval

    def __init__(self, function, interval):
        self._function = function
        self._interval = interval

    def evaluate(self, x: float) -> float:
        """ Evaluate the integrand """
        return self._function(x)

    def interval(self) -> IntegrandInterval:
        """ The integrand interval """
        return self._interval


__all__ = ["BaseIntegrand"]
