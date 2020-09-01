#!/usr/bin/env python3

from cProfile import Profile

from numpy import linspace

from super_material.gap_energy import BCSGapEnergy


def run():
    bcs_gap_energy = BCSGapEnergy(1.5e-3, 2.3)
    temperatures = linspace(0, bcs_gap_energy.critical_temperature(), 1000)

    with Profile() as profile:
        for temperature in temperatures:
            bcs_gap_energy.evaluate(temperature)

    profile.print_stats()


if __name__ == "__main__":
    run()
