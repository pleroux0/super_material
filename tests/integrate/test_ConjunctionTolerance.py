from material.integrate.ConjunctionTolerance import ConjunctionTolerance
from material.integrate.RelativeTolerance import RelativeTolerance
from material.integrate.AbsoluteTolerance import AbsoluteTolerance

from .test_ToleranceInterface import assert_tolerance_interface


def test_relative_tolerance():
    relative_tolerance = RelativeTolerance(1e-3)
    absolute_tolerance = AbsoluteTolerance(1e-3)

    tolerance_a = ConjunctionTolerance([absolute_tolerance])
    tolerance_r = ConjunctionTolerance([relative_tolerance])
    tolerance_ar = ConjunctionTolerance([absolute_tolerance, relative_tolerance])
    tolerance_ra = ConjunctionTolerance([relative_tolerance, absolute_tolerance])

    assert_tolerance_interface(tolerance_a)
    assert_tolerance_interface(tolerance_r)
    assert_tolerance_interface(tolerance_ar)
    assert_tolerance_interface(tolerance_ra)

    # All fail
    assert not tolerance_a.within_tolerance(1, 1.1)
    assert not tolerance_r.within_tolerance(1, 1.1)
    assert not tolerance_ar.within_tolerance(1, 1.1)
    assert not tolerance_ra.within_tolerance(1, 1.1)

    # Relative only
    assert not tolerance_a.within_tolerance(1000, 1000.5)
    assert tolerance_r.within_tolerance(1000, 1000.5)
    assert not tolerance_ar.within_tolerance(1000, 1000.5)
    assert not tolerance_ra.within_tolerance(1000, 1000.5)

    # Abs only
    assert tolerance_a.within_tolerance(1e-3, 1e-4)
    assert not tolerance_r.within_tolerance(1e-3, 1e-4)
    assert not tolerance_ar.within_tolerance(1e-3, 1e-4)
    assert not tolerance_ra.within_tolerance(1e-3, 1e-4)

    # Both
    assert tolerance_a.within_tolerance(1, 1.0005)
    assert tolerance_r.within_tolerance(1, 1.0005)
    assert tolerance_ar.within_tolerance(1, 1.0005)
    assert tolerance_ra.within_tolerance(1, 1.0005)
