from abc import ABC, abstractmethod

class SuperconductorConductivityInterface(ABC):

    @abstractmethod
    def evaluate(self, temperature: float, frequency: float) -> complex:
        """ Calculates the superconductor complex conductivity """
