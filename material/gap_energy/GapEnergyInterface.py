from abc import ABC, abstractmethod
from ..constants import h_bar, pi

class GapEnergyInterface(ABC):
    """ Interface to calculate the gap energy """

    @abstractmethod
    def evaluate(self, temperature: float) -> float:
        """ Evaluate the gap energy at the specific temperature """

    @abstractmethod
    def critical_temperature(self) -> float:
        """ Get the critical temperature of the gap energy """

    @abstractmethod
    def gap_energy_0(self) -> float:
        """ Get the gap energy at T = 0 K """

    def critical_frequency(self, temperature: float) -> float:
        return self.evaluate(temperature)/(h_bar * pi)

__all__ = ["GapEnergyInterface"]
