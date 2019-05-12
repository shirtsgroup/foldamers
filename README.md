This repository contains tools for building coarse-grained polymer models.  These models are designed to interface with OpenMM and PyRosetta.

## To use this repository install it using standard Python conventions:

python setup.py install

## Dependencies:

simtk.unit
simtk.openmm

## Getting started:

Test your installation by opening a new Python session and typing the following:

import foldamers

or

from foldamers import *

If this test does not work, please check to make sure that the foldamers path is included in the $PYTHONPATH system variable.  This is the most common place for the package installation to fail.

**For full documentation please refer to 'manual.pdf'.**
