from abc import ABC, abstractmethod


class SuperconductorConductivityInterface(ABC):
    """ Superconductor conductivity abstract class """

    @abstractmethod
    def evaluate(self, temperature: float, frequency: float) -> complex:
        """ Calculates the superconductor complex conductivity """


__all__ = ["SuperconductorConductivityInterface"]
