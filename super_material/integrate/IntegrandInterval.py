from .IntegrandBoundary import IntegrandBoundary


class IntegrandInterval:
    _start: IntegrandBoundary
    _end: IntegrandBoundary

    def __init__(self, start, end):
        self._start = start
        self._end = end

    def start(self) -> IntegrandBoundary:
        return self._start

    def end(self) -> IntegrandBoundary:
        return self._end


__all__ = ["IntegrandInterval"]
