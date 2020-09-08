# Super Material

Providing a constant interface to superconductor conductivity computations.

## Motivation

A lot of superconductor electromagnetic simulations can be transformed into classical computational electromagnetic (CEM) problems by substituting a complex conductivity term. Calculating these superconductor conductivity terms can be trivial, such as when calculating the two-fluid superconductor conductivity, but can also be numerically challenging, such as when calculating the Zimmermann superconductor conductivity. Calculating the superconductor conductivity term oneself does not give someone trying to solve a classical electromagnetic problem any advantage. With this library, one can directly jump into the CEM problem rather than spending time evaluating superconductor conductivity terms.

Other codes for calculating these terms are available, but some are numerically unstable and inefficient. We aim to provide a range of efficient and numerically stable set of routines.

## Installation

General users should install the latest release with pip

```bash
pip install super_material
```

Developers should install from the source directory using poetry

```bash
poetry install
```

## Usage example

A simple example showing how to calculate the Mattis-Bardeen superconductor conductivity for Niobium at 4.2 K and 100 GHz.

```python
from super_material import BCSGapEnergy, MattisBardeenSuperconductorConductivity

conductivity_0 = 2.4e7 # in Siemens per meter
temperature = 4.2 # in K
gap_energy_0 = 1.5e-3 # in eV
frequency = 100e9 # in Hz
kappa = 4000

gap_energy = BCSGapEnergy(gap_energy_0, 4000)
conductivity = MattisBardeenSuperconductorConductivity(gap_energy, conductivity_0)
result = conductivity.evaluate(temperature, frequency)
print(f"sigma = {result}")
```

For more information see the [full documentation](https://pleroux0.github.io/super_material/)

## Acknowledgements

This project was developed under IARPA contract SuperTools
(via the U.S. Army Research Office grant W911NF-17-1-0120).

## License

This project is licensed under the 2-Clause BSD license for maximum usability.

See [LICENSE.md](LICENSE.md) for more information
