from math import asin, sin, sqrt

from .IntegrandIntervalTransformInterface import IntegrandIntervalTransformInterface


class ChebyshevLowerSingularityTransform(IntegrandIntervalTransformInterface):
    _a: float

    def __init__(self, a: float):
        self._a = a

    def transform(self, x: float) -> float:
        return sqrt(x - self._a)

    def inverse_transform(self, u: float) -> float:
        return self._a + u ** 2

    def transform_derivative(self, x: float) -> float:
        return 0.5 / sqrt(x - self._a)


class ChebyshevUpperSingularityTransform(IntegrandIntervalTransformInterface):
    _b: float

    def __init__(self, b: float):
        self._b = b

    def transform(self, x: float) -> float:
        return sqrt(self._b - x)

    def inverse_transform(self, u: float) -> float:
        return self._b - u ** 2

    def transform_derivative(self, x: float) -> float:
        return -0.5 / sqrt(self._b - x)


class ChebyshevSingularityTransform(IntegrandIntervalTransformInterface):
    _a: float
    _b: float

    def __init__(self, a: float, b: float):
        self._a = a
        self._b = b

    def transform(self, x: float) -> float:
        return 2 * asin(sqrt((self._a - x) / (self._a - self._b)))

    def inverse_transform(self, u: float) -> float:
        return self._a + (self._b - self._a) * (sin(u / 2) ** 2)

    def transform_derivative(self, x: float) -> float:
        return 1 / sqrt((x - self._b) * (self._a - x))


__all__ = [
    "ChebyshevUpperSingularityTransform",
    "ChebyshevLowerSingularityTransform",
    "ChebyshevSingularityTransform",
]
