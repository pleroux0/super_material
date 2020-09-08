from abc import ABC, abstractmethod
from ..constants import h_bar, pi


class GapEnergyInterface(ABC):
    """ Superconductor gap energy abstract class """

    @abstractmethod
    def evaluate(self, temperature: float) -> float:
        """ Evaluate the gap energy at the specific temperature """

    @abstractmethod
    def critical_temperature(self) -> float:
        """ Get the critical temperature or transition temperature """

    @abstractmethod
    def gap_energy_0(self) -> float:
        """ Get the gap energy at T = 0 K """

    def critical_frequency(self, temperature: float) -> float:
        """ Get the critical frequency or gap frequency """
        return self.evaluate(temperature) / (h_bar * pi)


__all__ = ["GapEnergyInterface"]
