from abc import ABC, abstractmethod

from .IntegrandInterval import IntegrandInterval
from .IntegrandBoundary import IntegrandBoundary


class IntegrandIntervalTransformInterface(ABC):
    @abstractmethod
    def transform(self, x: float) -> float:
        """ Evaluates the transform """

    @abstractmethod
    def inverse_transform(self, u: float) -> float:
        """ Evaluates the inverse transform """

    @abstractmethod
    def transform_derivative(self, x: float) -> float:
        """ Evaluates the first derivative of the transform """

    def transform_boundary(self, boundary: IntegrandBoundary) -> IntegrandBoundary:
        transformed_value = self.transform(boundary.value())
        return IntegrandBoundary(transformed_value, boundary.defined_on_boundary())

    def transform_interval(self, interval: IntegrandInterval) -> IntegrandInterval:
        transformed_start = self.transform_boundary(interval.start())
        transformed_end = self.transform_boundary(interval.end())
        return IntegrandInterval(transformed_start, transformed_end)


__all__ = ["IntegrandIntervalTransformInterface"]
