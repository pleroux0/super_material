#!/usr/bin/env python3

import warnings
from math import pi

from numpy import array, linspace
import matplotlib.pyplot as plt

from super_material.conductivity import ZimmermannSuperconductorConductivity
from super_material.constants import h_bar
from super_material.gap_energy import BCSGapEnergy

warnings.simplefilter("error")


def calculate_normalized_zimmermann_conductivity(x: float, y: float):
    # Independant parameters
    conductivity_0 = 1
    gap_energy_0 = 1.5e-3
    kappa = 2.3
    temperature = 1e-6

    # Dependent parameters
    frequencies = x * gap_energy_0 / (h_bar * pi)
    scattering_time = h_bar / (2 * y * gap_energy_0)

    # Define gap energy and superconductor conductivity object
    gap_energy = BCSGapEnergy(gap_energy_0, kappa)
    conductivity = ZimmermannSuperconductorConductivity(
        gap_energy, conductivity_0, scattering_time
    )

    # Calculate conductivity values
    values = array([conductivity.evaluate(temperature, f) for f in frequencies])

    # Return the values
    return values


def zimmermann_example():
    # Determine frequency samples
    num_frequency_intervals = 7
    num_frequency_samples_per_interval = 20
    x = linspace(
        0,
        num_frequency_intervals,
        num_frequency_intervals * num_frequency_samples_per_interval
        + num_frequency_intervals
        + 1,
    )

    nozero_x = x[1:]

    # Define impurity parameters to plot
    impurity_parameters = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 500]
    num_impurity_parameters = len(impurity_parameters)

    # Calculate conductivity for each impurity parameter
    conductivity = [
        calculate_normalized_zimmermann_conductivity(nozero_x, y)
        for y in impurity_parameters
    ]

    # Setup figure
    _, (real_axis, imag_axis) = plt.subplots(ncols=1, nrows=2)

    # Setup x axis
    for axes in [real_axis, imag_axis]:
        axes.set_xlim(0, num_frequency_intervals)
        axes.set_xticks(range(0, num_frequency_intervals + 1))
        axes.set_xticklabels(
            ["0", "$\\omega_g$"]
            + [""] * (num_frequency_intervals - 2)
            + [f"${num_frequency_intervals} \\omega_g$"]
        )
        axes.set_xlabel("Frequency")

    # Setup y axis
    real_axis.set_ylabel("$\\Re\\{\\sigma_{sc}\\}$")
    imag_axis.set_ylabel("$\\Im\\{\\sigma_{sc}\\}$")

    for axes in [real_axis, imag_axis]:
        axes.set_ylim(0, 1)
        axes.set_yticks([0, 1])
        axes.set_yticklabels(["0", "$\\sigma_0$"])

    # Plot conductivity
    dashes = [
        (None, None),
        [2, 2, 2, 2],
        [4, 2, 4, 2],
        (None, None),
        [8, 2, 8, 2],
        [6, 2, 6, 2],
        [4, 2, 4, 2],
        [2, 2, 2, 2],
        (None, None),
        (None, None),
    ]

    double_dashes = [
        (None, None),
        [4, 4, 4, 4],
        [8, 4, 8, 4],
        (None, None),
        [16, 4, 16, 4],
        [12, 4, 12, 4],
        [8, 4, 8, 4],
        [4, 4, 4, 4],
        (None, None),
        (None, None),
    ]

    for index in range(num_impurity_parameters):
        s = conductivity[index]

        real_axis.plot(
            nozero_x[num_frequency_samples_per_interval:],
            s.real[num_frequency_samples_per_interval:],
            dashes=dashes[index],
            linestyle="--",
            color="black",
        )

        real_axis.plot(
            x,
            1 / (1 + x ** 2 / impurity_parameters[index] ** 2),
            dashes=double_dashes[index],
            linestyle="--",
            color="black",
            linewidth=0.5,
        )

        imag_axis.plot(
            nozero_x, s.imag, dashes=dashes[index], linestyle="--", color="black",
        )

    # Set figure spines
    for axes in [real_axis, imag_axis]:
        for spine in ["top", "right"]:
            axes.spines[spine].set_visible(False)

        for spine in ["left", "bottom"]:
            axes.spines[spine].set_color("gray")
            axes.spines[spine].set_linewidth(0.5)

    # Show the figure
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    zimmermann_example()
