======================================
Superconductor conductivity comparison
======================================

We look at the superconductor conductivity for niobium calculated with
different methods. In the dirty limit, the Zimmermann and Mattis-Bardeen
expression are almost indistinguishable. When including the scattering rate in
the Zimmermann superconductor conductivity, the graphs deviates from the dirty
case. From the comparison, we conclude that scattering has to be taken into
account for high-frequency problems.

.. image:: /static/images/SuperconductorConductivityComparison.svg
    :align: center

The code can be found in
``./examples/superconductor_conductivity_comparison.py`` and can be run with

.. code-block:: bash

    poetry run ./examples/superconductor_conductivity_comparison.py
