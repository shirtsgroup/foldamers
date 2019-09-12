Thermodynamic analysis tools for coarse grained modeling
========================================================

The 'foldamers' package contains wide-ranging tools to analyze simulation data in order to estimate the thermodynamic properties of a system.  In particular, the package leverages tools from `pymbar <https://pymbar.readthedocs.io/en/latest/index.html>`_ to enable analysis of and estimate expectation values for sampled and unsampled thermodynamic states.

Using MBAR to compute expectation values for structural & thermodynamic properties
----------------------------------------------------------------------------------

Shown below is the main 'foldamers' function used to re-weight simulation results with 'pymbar'.


.. automodule:: parameters.reweight
    :members: get_mbar_expectation

.. raw:: latex

    \newpage

Evaluating thermodynamic properties with 'pymbar'
--------------------------------------------------------

Shown below 'pymbar'-based tools/functions to evaluate thermodynamic properties within the 'foldamers' package.

.. automodule:: parameters.reweight
    :members: get_free_energy_differences, get_entropy_differences, get_enthalpy_differences

.. raw:: latex

    \newpage

Tools to evaluate thermodynamic properties with 'pymbar'
--------------------------------------------------------

Shown below are other 'pymbar'-based tools/functions that aid evaluation of thermodynamic properties within the 'foldamers' package.

.. automodule:: parameters.reweight
    :members: get_temperature_list, get_intermediate_temperatures, calc_temperature_spacing, get_decorrelated_samples

.. raw:: latex

    \newpage


Calculating the heat capacity with pymbar
-----------------------------------------

Shown below is the primary function used to evaluate the heat capacity with pymbar:

.. automodule:: thermo.calc
    :members: get_heat_capacity

Shown below are functions/tools that can be used in order to calculate the heat capacity with pymbar.

.. automodule:: thermo.calc
    :members:
    :exclude-members: get_heat_capacity

