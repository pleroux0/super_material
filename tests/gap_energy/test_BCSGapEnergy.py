from dataclasses import dataclass

import numpy as np

from super_material.gap_energy.BCSGapEnergy import BCSGapEnergy

from .test_GapEnergyInterface import assert_gap_energy_interface


@dataclass
class BCSGapEnergyTestCase:
    name: str
    gap_energy_0: float
    kappa: float
    temperatures: np.ndarray
    expected: np.ndarray


def assert_bcs_gap_energy_test_case(test_case: BCSGapEnergyTestCase):
    gap_energy = BCSGapEnergy(test_case.gap_energy_0, test_case.kappa)

    # Interface
    assert_gap_energy_interface(gap_energy)

    # Expected evaluation
    data = [gap_energy.evaluate(temperature) for temperature in test_case.temperatures]
    assert np.allclose(data, test_case.expected)


def test_bcs_gap_energy():
    niobium_test_case = BCSGapEnergyTestCase(
        "Niobium",
        1.5e-3,
        2.3,
        np.array([2.1, 4.2, 6.3, 8.4]),
        np.array(
            [
                0.00149965137842116,
                0.001466898544413353,
                0.0013037620136030783,
                0.0008657797279237159,
            ]
        ),
    )

    assert_bcs_gap_energy_test_case(niobium_test_case)
