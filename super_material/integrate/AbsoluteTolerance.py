from .ToleranceInterface import ToleranceInterface


class AbsoluteTolerance(ToleranceInterface):
    _absolute_tolerance: float

    def __init__(self, absolute_tolerance: float):
        self._absolute_tolerance = absolute_tolerance

    def absolute_tolerance(self) -> float:
        return self._absolute_tolerance

    def within_tolerance(self, value: float, reference: float) -> bool:
        absolute_difference = abs(reference - value)
        return absolute_difference < self._absolute_tolerance


__all__ = ["AbsoluteTolerance"]
