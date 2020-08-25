from material.integrate import ToleranceInterface


def assert_tolerance_interface(tolerance: ToleranceInterface):
    # Should always be True
    assert tolerance.within_tolerance(1, 1)
