Utilities for the 'foldamers' package
=====================================

The 'foldamers' package contains a number of utilities for reading and writing data, as well as random structure generation.

Random structure generation in 'foldamers'
------------------------------------------

Random structures are often needed as a starting point for simulation of a new coarse grained model.  Shown below are the main tools that allows the user to build a random structure:

.. automodule:: utilities.util
    :members: random_positions, get_structure_from_library

.. raw:: latex

    \newpage

Shown below are other tools that support the task of random structure generation:

.. automodule:: utilities.util
    :members:
    :exclude-members: random_positions, get_structure_from_library

.. raw:: latex

    \newpage

Reading and writing to output within the 'foldamers' package
--------------------------------------------------------

Shown below are 'foldamers' tools for writing data and structures to output files.

.. automodule:: utilities.iotools
    :members:
