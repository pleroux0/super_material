from abc import ABC, abstractmethod

from typing import List


class FixedQuadratureInterface(ABC):
    @abstractmethod
    def weights(self) -> List[float]:
        """ Returns the weights of the quadratures points """

    @abstractmethod
    def abscissae(self) -> List[float]:
        """ Returns the weights of the quadrature points """

    @abstractmethod
    def num_quadrature_points(self) -> int:
        """ The number of quadrature points in the fixed quadrature """
