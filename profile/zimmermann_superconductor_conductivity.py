#!/usr/bin/env python3

from cProfile import Profile

from numpy import linspace

from super_material.gap_energy import BCSGapEnergy
from super_material.conductivity import ZimmermannSuperconductorConductivity


def run():
    gap_energy = BCSGapEnergy(1.5e-3, 4000)
    conductivity = ZimmermannSuperconductorConductivity(gap_energy, 2.4e7, 3e-14)

    frequencies = linspace(10e9, 1500e9, 200)

    with Profile() as profile:
        for frequency in frequencies:
            conductivity.evaluate(4.2, frequency)

    profile.print_stats()


if __name__ == "__main__":
    run()
