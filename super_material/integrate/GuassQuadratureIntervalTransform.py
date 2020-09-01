from .LinearIntegrandIntervalTransform import LinearIntegrandIntervalTransform
from .IntegrandInterval import IntegrandInterval


class LinearGuassQuadratureIntervalTransform(LinearIntegrandIntervalTransform):
    _base: IntegrandInterval

    def __init__(self, base: IntegrandInterval):
        assert base.start().is_finite()
        assert base.end().is_finite()

        start = base.start().value()
        end = base.end().value()

        m = 2 / (end - start)
        c = (start + end) / (start - end)

        super().__init__(m, c)

        self._base = base

    def base(self) -> IntegrandInterval:
        return self._base


__all__ = ["LinearGuassQuadratureIntervalTransform"]
