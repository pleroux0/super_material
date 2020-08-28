#!/usr/bin/env python3

from numpy import linspace, array
from matplotlib.pyplot import subplots, show

from material.gap_energy import BCSGapEnergy


def bcs_gap_energy_example():
    # Niobium paarameters
    gap_energy_0 = 1.5e-3
    kappa = 2.3

    # Define gap energy object
    bcs_gap_energy = BCSGapEnergy(gap_energy_0, kappa)

    # Determine temperatures
    Tc = bcs_gap_energy.critical_temperature()
    temperatures = linspace(0, Tc, 200)

    # Calculate gap energy at specified temperatures
    values = array([bcs_gap_energy.evaluate(T) for T in temperatures])

    # Setup figure
    _, axes = subplots()

    # Plot normalized gap energy vs temperature
    axes.plot(temperatures / Tc, values / gap_energy_0, color="black")

    # Setup x axis
    axes.set_xlim(0, 1 + 1e-3)
    axes.set_xticks([0, 1])
    axes.set_xticklabels(["0", "$T_c$"])
    axes.set_xlabel("Temperature")

    # Setup y axis
    axes.set_ylim(0, 1 + 5e-3)
    axes.set_yticks([0, 1])
    axes.set_yticklabels(["0", "$\\Delta_0$"])
    axes.set_ylabel("Gap energy")

    # Set figure spines
    for spine in ["top", "right"]:
        axes.spines[spine].set_visible(False)

    for spine in ["left", "bottom"]:
        axes.spines[spine].set_color("gray")
        axes.spines[spine].set_linewidth(0.5)

    # Show the figure
    show()


if __name__ == "__main__":
    bcs_gap_energy_example()
