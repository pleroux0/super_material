from abc import ABC, abstractmethod

from .IntegrandInterval import IntegrandInterval


class IntegrandInterface(ABC):
    @abstractmethod
    def evaluate(self, x: float) -> float:
        """ Evaluate the integrand """

    @abstractmethod
    def interval(self) -> IntegrandInterval:
        """ Evaluate the integrand """


__all__ = ["IntegrandInterface"]
