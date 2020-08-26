from math import isclose

from material.gap_energy.GapEnergyInterface import GapEnergyInterface


def assert_gap_energy_interface(gap_energy: GapEnergyInterface):
    Tc = gap_energy.critical_temperature()
    Delta0 = gap_energy.gap_energy_0()

    assert isclose(gap_energy.evaluate(Tc), 0)
    assert isclose(gap_energy.evaluate(0), Delta0)
