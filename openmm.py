#!/usr/local/bin/env python

# ==============================================================================
# MODULE DOCSTRING
# ==============================================================================

"""
foldamers/openmm.py
====

Tools for building coarse grained models 
and running simulations with OpenMM/Yank.

"""

# ==============================================================================
# GLOBAL IMPORTS
# ==============================================================================


from simtk import openmm as mm
from simtk import unit
import simtk.openmm.app.element as elem

def get_box_vectors(box_size):
        """

        Construct an OpenMM system object for our coarse grained model

        Parameters
        ----------

        box_size: Simulation box length ( float * simtk.unit.length )

        """

        units = box_size.unit
        a = unit.Quantity(np.zeros([3]), units)
        a[0] = box_size
        b = unit.Quantity(np.zeros([3]), units)
        b[1] = box_size
        c = unit.Quantity(np.zeros([3]), units)
        c[2] = box_size
        return([a,b,c])

def set_box_vectors(system,box_size):
        """

        Construct an OpenMM system object for our coarse grained model

        Parameters
        ----------

        system: OpenMM system object

        box_size: Simulation box length ( float * simtk.unit.length )

        """

        a,b,c = get_box_vectors(box_size)
        system.setDefaultPeriodicBoxVectors(a, b, c)
        return(system)

def build_mm_system(box_size,mass,num_beads):
        """

        Construct an OpenMM system object for our coarse grained model

        Parameters
        ----------

        box_size: Simulation box length ( float * simtk.unit.length )

        mass: Coarse grained particle mass ( float * simtk.unit.length )

        num_beads: Total number of beads in our coarse grained model (int)

        """

        system = mm.System()

        for particle in range(num_beads):
          system.addParticle(mass)

        return(system)

def build_mm_topology(num_beads):
        """

        Construct an OpenMM system object for our coarse grained model

        Parameters
        ----------

        box_size: Simulation box length ( float * simtk.unit.length )

        mass: Coarse grained particle mass ( float * simtk.unit.length )

        num_beads: Total number of beads in our coarse grained model (int)

        """

        topology = mm.Topology()

        for particle in range(num_beads):
          system.addParticle(mass)

        return(topology)

class cgmodel(object):
        """
        Construct all of the objects that OpenMM expects/requires 
        for simulations with a coarse grained model.

        Parameters
        ----------

        box_size: Simulation box length, 
        default = 10.00 * unit.nanometer

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

        box_size
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
        topology
        system

        Notes
        -----
        
        """

        # Built in class attributes
        _BUILT_IN_REGIONS = ('box_size','polymer_length','backbone_length','sidechain_length','sidechain_positions','mass','sigma','epsilon','bond_length','bs_bond_length','bb_bond_length','ss_bond_length','charge','topology','system')

        def __init__(self, box_size = 10.00 * unit.nanometer, polymer_length = 12, backbone_length = 1, sidechain_length = 1, sidechain_positions = [0], mass = 12.0 * unit.amu, sigma = 8.4 * unit.angstrom, epsilon = 0.5 * unit.kilocalorie_per_mole, bond_length = 1.0 * unit.angstrom, bb_bond_length = 1.0 * unit.angstrom, bs_bond_length = 1.0 * unit.angstrom, ss_bond_length = 1.0 * unit.angstrom, charge = 0.0 * unit.elementary_charge):

          """
          Initialize variables that were passed as input
          """

          self._box_size = box_size
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

          self._num_beads = polymer_length * ( backbone_length + sidechain_length )
          self._system = build_mm_system(box_size,mass,self._num_beads)

        """
        Initialize attributes of our coarse grained model.
        """

        @property
        def box_size(self):
          """Make the 'box_size' a property of this 'cgmodel' object."""
          return self._box_size

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
        def topology(self)
          return self._topology

        @property
        def system(self)
          return self._system
          
