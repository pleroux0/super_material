from typing import List

from .ToleranceInterface import ToleranceInterface


class DisjunctionTolerance(ToleranceInterface):
    _children: List[ToleranceInterface]

    def __init__(self, children: List[ToleranceInterface]):
        self._children = children

    def children(self) -> List[ToleranceInterface]:
        return self._children

    def within_tolerance(self, value: float, reference: float) -> bool:
        for child in self._children:
            if child.within_tolerance(value, reference):
                return True

        return False


__all__ = ["DisjunctionTolerance"]
