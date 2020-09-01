from math import isfinite, isnan


class IntegrandBoundary:
    _value: float
    _defined_on_boundary: bool

    def __init__(self, value, defined_on_boundary: bool):
        assert not isnan(value)

        self._value = value
        self._defined_on_boundary = defined_on_boundary

    def is_finite(self) -> bool:
        return isfinite(self.value())

    def value(self) -> float:
        return self._value

    def defined_on_boundary(self) -> bool:
        return self._defined_on_boundary


__all__ = ["IntegrandBoundary"]
