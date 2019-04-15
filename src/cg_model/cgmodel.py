#!/usr/local/bin/env python

#Tools for building coarse grained models 

# ==============================================================================
# GLOBAL IMPORTS
# ==============================================================================

from simtk import unit
import sys, os
sys.path.insert(0, os.path.abspath('/home/gmeek/Foldamers/foldamers/src'))
from utilities import util as util

class CGModel(object):
        """
        Construct a coarse grained model.

        Parameters
        ----------

        polymer_length: Number of monomer units (integer), default = 8
      
        backbone_length: Number of beads in the backbone 
        portion of each (individual) monomer (integer), default = 1

        sidechain_length: Number of beads in the sidechain
        portion of each (individual) monomer (integer), default = 1

        sidechain_positions: List of integers defining the backbone
        bead indices upon which we will place the sidechains,
        default = [0] (Place a sidechain on the backbone bead with
        index "0" (first backbone bead) in each (individual) monomer

        mass: Mass of coarse grained beads ( float * simtk.unit.mass )
        default = 12.0 * unit.amu

        sigma: Non-bonded bead Lennard-Jones interaction distances,
        ( float * simtk.unit.distance )
        default = 8.4 * unit.angstrom

        epsilon: Non-bonded bead Lennard-Jones interaction strength,
        ( float * simtk.unit.energy )
        default = 0.5 * unit.kilocalorie_per_mole

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        bb_bond_length: Bond length for all bonded backbone beads,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        bs_bond_length: Bond length for all backbone-sidechain bonds,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        ss_bond_length: Bond length for all beads within a sidechain,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        charge: Charge for all beads
        ( float * simtk.unit.charge )
        default = 0.0 * unit.elementary_charge

        Attributes
        ----------

        polymer_length
        backbone_length
        sidechain_length
        sidechain_positions
        mass
        sigma
        epsilon
        bond_length
        bb_bond_length
        bs_bond_length
        ss_bond_length
        charge
        num_beads
        positions
        topology

        Notes
        -----
        
        """

        # Built in class attributes
        _BUILT_IN_REGIONS = ('polymer_length','backbone_length','sidechain_length','sidechain_positions','mass','sigma','epsilon','bond_length','bs_bond_length','bb_bond_length','ss_bond_length','charge','num_beads','positions','topology')

        def __init__(self, polymer_length = 12, backbone_length = 1, sidechain_length = 1, sidechain_positions = [0], mass = 12.0 * unit.amu, sigma = 8.4 * unit.angstrom, epsilon = 0.5 * unit.kilocalorie_per_mole, bond_length = 1.0 * unit.angstrom, bb_bond_length = 1.0 * unit.angstrom, bs_bond_length = 1.0 * unit.angstrom, ss_bond_length = 1.0 * unit.angstrom, charge = 0.0 * unit.elementary_charge):

          """
          Initialize variables that were passed as input
          """

          self._polymer_length = polymer_length
          self._backbone_length = backbone_length
          self._sidechain_length = sidechain_length
          self._sidechain_positions = sidechain_positions
          self._mass = mass
          self._sigma = sigma
          self._epsilon = epsilon
          self._bond_length = bond_length
          self._bb_bond_length = bb_bond_length
          self._bs_bond_length = bs_bond_length
          self._ss_bond_length = ss_bond_length
          self._charge = charge         

          """
          Initialize new (coarse grained) particle types:
          """

          self._num_beads = polymer_length * ( backbone_length + sidechain_length )

          self._positions = util.random_positions( polymer_length, backbone_length, sidechain_length, sidechain_positions, bond_length, sigma ) 

          """
          Initialize attributes of our coarse grained model.
          """

        @property
        def polymer_length(self):
          return self._polymer_length

        @property
        def backbone_length(self):
          return self._backbone_length

        @property
        def sidechain_length(self):
          return self._sidechain_length

        @property
        def sidechain_positions(self):
          return self._sidechain_positions

        @property
        def mass(self):
          return self._mass

        @property
        def sigma(self):
          return self._sigma

        @property
        def epsilon(self):
          return self._epsilon

        @property
        def bond_length(self):
          return self._bond_length

        @property
        def bb_bond_length(self):
          return self._bb_bond_length

        @property
        def bs_bond_length(self):
          return self._bs_bond_length

        @property
        def ss_bond_length(self):
          return self._ss_bond_length

        @property
        def charge(self):
          return self._charge

        @property
        def num_beads(self):
          return self._num_beads

        @property
        def topology(self):
          return self._topology

