from dataclasses import dataclass

import numpy as np

from material.conductivity.ZimmermannSuperconductorConductivity import (
    ZimmermannSuperconductorConductivity,
)

from material.gap_energy import BCSGapEnergy


@dataclass
class ZimmermannSuperconductorConductivityTestCase:
    name: str
    gap_energy_0: float
    kappa: float
    temperature: float
    scattering_time: float
    conductivity_0: float
    frequencies: np.ndarray
    expected: np.ndarray


def assert_zimmermann_superconductor_conductivity_test_case(
    test_case: ZimmermannSuperconductorConductivityTestCase,
):
    gap_energy = BCSGapEnergy(test_case.gap_energy_0, test_case.kappa)
    conductivity = ZimmermannSuperconductorConductivity(
        gap_energy, test_case.conductivity_0, test_case.scattering_time
    )

    # Expected evaluation
    data = [
        conductivity.evaluate(test_case.temperature, frequency)
        for frequency in test_case.frequencies
    ]
    print(data)
    assert np.allclose(data, test_case.expected)


def test_zimmermann_superconductor_conductivity():
    niobium_4_2K_test_case = ZimmermannSuperconductorConductivityTestCase(
        "Niobium 4.2K",
        1.5e-3,
        2.3,
        4.2,
        3e-14,
        2.4e7,
        [10e9, 50e9, 100e9, 200e9, 400e9, 700e9, 750e9, 900e9, 1500e9],
        [
            (10250884.753168825 + 1972166271.6813483j),
            (5161389.872991417 + 397268365.18754834j),
            (3219907.764522122 + 199399731.12885547j),
            (1692024.181196706 + 99127792.46029282j),
            (753010.5041757216 + 46736249.72993765j),
            (371955.83095510636 + 18918296.25965787j),
            (2295643.7092327923 + 14733203.022071129j),
            (7630102.486317393 + 9528419.61700693j),
            (15999891.283731477 + 6190404.8912611175j),
        ],
    )

    assert_zimmermann_superconductor_conductivity_test_case(niobium_4_2K_test_case)
