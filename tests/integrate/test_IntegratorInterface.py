from material.integrate import IntegrandInterface, IntegrandBoundary, IntegrandInterval


class ParabolicTestIntegrand(IntegrandInterface):
    def evaluate(self, x: float) -> float:
        return x ** 2

    def interval(self) -> IntegrandInterval:
        start = IntegrandBoundary(0, True)
        end = IntegrandBoundary(4, True)
        interval = IntegrandInterval(start, end)
        return interval

    @staticmethod
    def analytical() -> float:
        return 4 ** 3 / 3
