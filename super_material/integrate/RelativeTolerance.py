from .ToleranceInterface import ToleranceInterface


class RelativeTolerance(ToleranceInterface):
    _relative_tolerance: float

    def __init__(self, relative_tolerance: float):
        self._relative_tolerance = relative_tolerance

    def relative_tolerance(self) -> float:
        return self._relative_tolerance

    def within_tolerance(self, value: float, reference: float) -> bool:
        absolute_difference = abs(reference - value)
        relative_difference = abs(absolute_difference / reference)
        return relative_difference < self._relative_tolerance


__all__ = ["RelativeTolerance"]
