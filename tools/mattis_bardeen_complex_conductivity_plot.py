#!env python

import warnings
from math import pi

import numpy as np
import matplotlib.pyplot as plt

from material.conductivity import MattisBardeenComplexConductivity
from material.constants import h_bar
from material.gap_energy import BCSGapEnergy

warnings.simplefilter("error")

def miniapp(
    conductivity_0: float, gap_energy_0: float, kappa: float, temperature: float
):
    gap_energy = BCSGapEnergy(gap_energy_0, kappa)
    gap_energy_T = gap_energy.evaluate(temperature)
    conductivity = MattisBardeenComplexConductivity(gap_energy, conductivity_0)

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

    plt.figure()
    plt.plot(frequency_values, np.abs(conductivity_values.real))
    plt.plot(frequency_values, np.abs(conductivity_values.imag))
    plt.xlim(0, max_frequency)
    plt.ylim(0, conductivity_0)
    plt.show()


if __name__ == "__main__":
    miniapp(2.4e7, 1.5e-3, 2.3, 4.2)
