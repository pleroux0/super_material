from .IntegrandIntervalTransformInterface import IntegrandIntervalTransformInterface


class LinearIntegrandIntervalTransform(IntegrandIntervalTransformInterface):
    _m: float
    _c: float

    def __init__(self, m: float, c: float):
        self._m = m
        self._c = c

    def m(self):
        return self._m

    def c(self):
        return self._c

    def transform(self, x: float) -> float:
        return self.m() * x + self.c()

    def inverse_transform(self, u: float) -> float:
        return (u - self.c()) / self.m()

    def transform_derivative(self, x: float) -> float:
        return self.m()


__all__ = ["LinearIntegrandIntervalTransform"]
