from typing import List

from .ToleranceInterface import ToleranceInterface


class ConjunctionTolerance(ToleranceInterface):
    _children: List[ToleranceInterface]

    def __init__(self, children: List[ToleranceInterface]):
        self._children = children

    def children(self) -> List[ToleranceInterface]:
        return self._children

    def within_tolerance(self, value: float, reference: float) -> bool:
        for child in self._children:
            if not child.within_tolerance(value, reference):
                return False

        return True


__all__ = ["ConjunctionTolerance"]
