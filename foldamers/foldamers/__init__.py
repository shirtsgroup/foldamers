"""
Foldamers
This repository contains tools for building coarse-grained polymer models. These models are designed to interface with OpenMM and PyRosetta.
"""

# Add imports here
from .foldamers import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions
