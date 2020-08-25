from abc import ABC, abstractmethod

class SuperconductorComplexConductivityInterface(ABC):

    @abstractmethod
    def evaluate(self, temperature: float, frequency: float):
        """ Calculates the superconductor complex conductivity """
