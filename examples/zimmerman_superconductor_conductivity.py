#!/usr/bin/env python3

import warnings
from math import pi

import numpy as np
import matplotlib.pyplot as plt

from material.conductivity import ZimmermannSuperconductorConductivity
from material.constants import h_bar
from material.gap_energy import BCSGapEnergy

warnings.simplefilter("error")


def miniapp(
        conductivity_0: float, gap_energy_0: float, kappa: float, temperature: float, scattering_time: float
):
    gap_energy = BCSGapEnergy(gap_energy_0, kappa)
    gap_energy_T = gap_energy.evaluate(temperature)
    conductivity = ZimmermannSuperconductorConductivity(gap_energy, conductivity_0, scattering_time)

    gap_frequency = gap_energy_T / (pi * h_bar)
    n = 4
    max_frequency = gap_frequency * n

    frequency_values = np.linspace(0, max_frequency, (n + 2) + (n + 1) * 100)[1:]
    conductivity_values = np.array(
        [
            conductivity.evaluate(temperature, frequency)
            for frequency in frequency_values
        ]
    )

    frequency_scale = gap_frequency
    conductivity_scale = conductivity_0

    plt.figure()

    plt.plot(
        frequency_values / frequency_scale,
        conductivity_values.real / conductivity_scale,
    )

    plt.plot(
        frequency_values / frequency_scale,
        conductivity_values.imag / conductivity_scale,
    )

    plt.xlabel("Frequency [$\\omega_g$]")
    plt.xlim(0, max_frequency / frequency_scale)

    plt.ylabel("Conductivity [$\\sigma_0$]")
    plt.ylim(0, conductivity_0 / conductivity_scale)

    plt.show()


if __name__ == "__main__":
    miniapp(2.4e7, 1.5e-3, 2.3, 4.2, 3e-14*1e-6)
