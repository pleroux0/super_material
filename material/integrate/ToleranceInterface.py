""" Integration error tolerance """

from abc import ABC, abstractmethod


class ToleranceInterface(ABC):
    """ Interface to define integration error tolerances """

    @abstractmethod
    def within_tolerance(self, value: float, reference: float) -> bool:
        """ Returns wether the value is within the tolerance from the reference """


__all__ = ["ToleranceInterface"]
