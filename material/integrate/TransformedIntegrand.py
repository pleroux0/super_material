from .IntegrandInterface import IntegrandInterface
from .IntegrandIntervalTransformInterface import IntegrandIntervalTransformInterface


class TransformedIntegrand(IntegrandInterface):
    _base: IntegrandInterface
    _transform: IntegrandIntervalTransformInterface

    def __init__(
        self, base: IntegrandInterface, transform: IntegrandIntervalTransformInterface
    ):
        self._base = base
        self._transform = transform

    def evaluate(self, u: float) -> float:
        x = self._transform.inverse_transform(u)
        scale = self._transform.transform_derivative(x)

        return self._base.evaluate(x) / scale

    def interval(self) -> IntegrandInterface:
        base_interval = self._base.interval()
        return self._transform.transform_interval(base_interval)


__all__ = ["TransformedIntegrand"]
