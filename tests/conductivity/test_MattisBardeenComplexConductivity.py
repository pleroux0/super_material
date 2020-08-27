from dataclasses import dataclass

import numpy as np

from material.conductivity.MattisBardeenComplexConductivity import (
    MattisBardeenComplexConductivity,
)

from material.gap_energy import BCSGapEnergy


@dataclass
class MattisBardeenSuperconductorConductivityTestCase:
    name: str
    gap_energy_0: float
    kappa: float
    temperature: float
    conductivity_0: float
    frequencies: np.ndarray
    expected: np.ndarray


def assert_mattis_bardeen_superconductor_conductivity_test_case(
    test_case: MattisBardeenSuperconductorConductivityTestCase,
):
    gap_energy = BCSGapEnergy(test_case.gap_energy_0, test_case.kappa)
    conductivity = MattisBardeenComplexConductivity(
        gap_energy, test_case.conductivity_0
    )

    # Expected evaluation
    data = [
        conductivity.evaluate(test_case.temperature, frequency)
        for frequency in test_case.frequencies
    ]
    assert np.allclose(data, test_case.expected)


def test_mattis_bardeen_complex_conductivity():
    niobium_4_2K_test_case = MattisBardeenSuperconductorConductivityTestCase(
        "Niobium 4.2K",
        1.5e-3,
        2.3,
        4.2,
        2.4e7,
        [10e9, 50e9, 100e9, 200e9, 400e9, 700e9, 750e9, 900e9, 1500e9],
        [
            (10265228.751330435 + 2587824292.7723484j),
            (5175645.913171842 + 520232687.85833454j),
            (3233698.424668898 + 260619985.0545055j),
            (1705037.597273664 + 129211990.12798274j),
            (765454.21826278 + 60713758.56712125j),
            (384935.2244671815 + 25161604.419847064j),
            (2311924.952897775 + 20180741.739129208j),
            (7707437.14390103 + 12955217.875369202j),
            (16873019.942420434 + 4344969.773031883j),
        ],
    )

    assert_mattis_bardeen_superconductor_conductivity_test_case(niobium_4_2K_test_case)
