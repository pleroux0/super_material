#!/usr/bin/env python3

import warnings
from math import pi

import numpy as np
import matplotlib.pyplot as plt

from super_material.conductivity import (
    MattisBardeenSuperconductorConductivity,
    ZimmermannSuperconductorConductivity,
)
from super_material.constants import h_bar
from super_material.gap_energy import BCSGapEnergy

warnings.simplefilter("error")


def superconductor_conductivity_comparison_example():
    # Niobium paarameters
    gap_energy_0 = 1.5e-3
    kappa = 2.3
    conductivity_0 = 2.4e7
    scattering_time = 3e-14
    dirty_scattering_time = 0
    temperature = 4.2

    # Create gap energy object
    gap_energy = BCSGapEnergy(gap_energy_0, kappa)

    # Create superconductor conductivity objects
    mattis_bardeen = MattisBardeenSuperconductorConductivity(gap_energy, conductivity_0)

    zimmerman_dirty = ZimmermannSuperconductorConductivity(
        gap_energy, conductivity_0, dirty_scattering_time
    )

    zimmerman = ZimmermannSuperconductorConductivity(
        gap_energy, conductivity_0, scattering_time
    )

    # Setup frequency information
    num_frequencies = 501
    maximum_frequency = 1500e9
    frequencies = np.linspace(0, maximum_frequency, num_frequencies + 1)[1:]

    # Calculate conductivity over frequency
    conductivities = []

    for interface in [mattis_bardeen, zimmerman_dirty, zimmerman]:
        s = np.array([interface.evaluate(temperature, f) for f in frequencies])
        conductivities.append(s)

    # Setup figure
    _, (real_axis, imag_axis) = plt.subplots(ncols=1, nrows=2)

    # Setup x axis
    for axes in [real_axis, imag_axis]:
        axes.set_xlim(0, maximum_frequency)
        axes.set_xlabel("Frequency")

    # Setup y axis
    real_axis.set_ylabel("$\\Re\\{\\sigma_{sc}\\}$")
    imag_axis.set_ylabel("$\\Im\\{\\sigma_{sc}\\}$")

    for axes in [real_axis, imag_axis]:
        axes.set_ylim(0, conductivity_0)
        axes.set_yticks([0, conductivity_0])
        axes.set_yticklabels(["0", "$\\sigma_0$"])

    # Plot conductivities
    dashes = [
        (None, None),
        [4, 2, 4, 2],
        [8, 2, 1, 2],
    ]

    for index, s in enumerate(conductivities):
        real_axis.plot(
            frequencies, s.real, linestyle="--", color="black", dashes=dashes[index]
        )

        imag_axis.plot(
            frequencies, s.imag, linestyle="--", color="black", dashes=dashes[index]
        )

    # Set figure spines
    for axes in [real_axis, imag_axis]:
        for spine in ["top", "right"]:
            axes.spines[spine].set_visible(False)

        for spine in ["left", "bottom"]:
            axes.spines[spine].set_color("gray")
            axes.spines[spine].set_linewidth(0.5)

    # Show figure
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    superconductor_conductivity_comparison_example()
