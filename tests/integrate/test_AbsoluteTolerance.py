from material.integrate.AbsoluteTolerance import AbsoluteTolerance

from .test_ToleranceInterface import assert_tolerance_interface


def test_relative_tolerance():
    tolerance_1 = AbsoluteTolerance(1e-1)
    tolerance_2 = AbsoluteTolerance(1e-2)
    tolerance_3 = AbsoluteTolerance(1e-3)

    assert_tolerance_interface(tolerance_1)
    assert_tolerance_interface(tolerance_2)
    assert_tolerance_interface(tolerance_3)

    assert not tolerance_1.within_tolerance(1, 1.5)
    assert not tolerance_2.within_tolerance(1, 1.5)
    assert not tolerance_3.within_tolerance(1, 1.5)

    assert tolerance_1.within_tolerance(1, 1.05)
    assert not tolerance_2.within_tolerance(1, 1.05)
    assert not tolerance_3.within_tolerance(1, 1.05)

    assert tolerance_1.within_tolerance(1, 1.005)
    assert tolerance_2.within_tolerance(1, 1.005)
    assert not tolerance_3.within_tolerance(1, 1.005)

    assert tolerance_1.within_tolerance(1, 1.0005)
    assert tolerance_2.within_tolerance(1, 1.0005)
    assert tolerance_3.within_tolerance(1, 1.0005)
