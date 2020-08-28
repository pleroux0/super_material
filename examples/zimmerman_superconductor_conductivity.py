#!/usr/bin/env python3

import warnings
from math import pi

import numpy as np
import matplotlib.pyplot as plt

from material.conductivity import ZimmermannSuperconductorConductivity
from material.constants import h_bar
from material.gap_energy import BCSGapEnergy

warnings.simplefilter("error")


def zimmermann_example():
    conductivity_0 = 1
    gap_energy_0 = 1.5e-3
    kappa = 2.3
    temperature = 1e-3
    gap_energy = BCSGapEnergy(gap_energy_0, kappa)

    n = 50
    x = np.linspace(0, 7, 7 * n + 6)[1:]
    omega = (2 * x * gap_energy_0) / h_bar
    frequency = omega / (2 * np.pi)

    y = np.array([0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 500])

    conductivity = []

    for yi in y:
        scattering_time = h_bar / (2 * yi * gap_energy_0)
        zimmermann = ZimmermannSuperconductorConductivity(
            gap_energy, conductivity_0, scattering_time
        )
        s = [zimmermann.evaluate(temperature, f) for f in frequency]
        conductivity.append(np.array(s))

    _, axes = plt.subplots(nrows=2, ncols=1)

    for index in range(len(y)):
        s = conductivity[index]
        axes[0].plot(x, s.real)
        axes[1].plot(x, s.imag)

    axes[0].set_xlim(0, 7)
    axes[1].set_xlim(0, 7)
    axes[0].set_ylim(0, 1)
    axes[1].set_ylim(0, 1)

    plt.show()


if __name__ == "__main__":
    zimmermann_example()
