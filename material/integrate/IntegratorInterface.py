from abc import ABC, abstractmethod

from .IntegrandInterface import IntegrandInterface


class IntegratorInterface(ABC):
    @abstractmethod
    def integrate(self, integrand: IntegrandInterface) -> float:
        """ Evaluate the definite integral """


__all__ = ["IntegratorInterface"]
