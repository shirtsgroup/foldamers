#!/usr/local/bin/env python

#Tools for building coarse grained models 

# ==============================================================================
# GLOBAL IMPORTS
# ==============================================================================

from simtk import unit
import sys, os
from ..utilities import util
from simtk import openmm as mm
import simtk.openmm.app.element as elem

def get_particle_masses(cgmodel):
        list_of_masses = []
        for backbone_bead in range(cgmodel.backbone_length):
            list_of_masses.append(cgmodel.mass)
            if backbone_bead in cgmodel.sidechain_positions:
              for sidechain in range(cgmodel.sidechain_length):
                 list_of_masses.append(cgmodel.mass)
        return(list_of_masses)

def add_new_elements(cgmodel,list_of_masses):
        element_index = 117
        mass_index = 0
        cg_particle_index = 1
        for backbone_bead in range(cgmodel.backbone_length):

         particle_name = str("bb-"+str(cg_particle_index))
         particle_symbol = str("B"+str(cg_particle_index))
         if particle_symbol not in elem.Element._elements_by_symbol:
          elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
          element_index = element_index + 1
          cg_particle_index = cg_particle_index + 1
          mass_index = mass_index + 1
         if backbone_bead in cgmodel.sidechain_positions:
           for sidechain in range(cgmodel.sidechain_length):
            if particle_symbol not in elem.Element._elements_by_symbol:
             particle_name = str("sc-"+str(cg_particle_index))
             particle_symbol = str("S"+str(cg_particle_index))
             elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
             element_index = element_index + 1
             cg_particle_index = cg_particle_index + 1
             mass_index = mass_index + 1

        return

def get_parent_bead(cgmodel,bead_index,backbone_bead_index=None,sidechain_bead=False):
              parent_bead = -1
              if bead_index != 1:
               if sidechain_bead == True:
                parent_bead = bead_index - 1
                return(parent_bead)
               else:
                 if backbone_bead_index - 2 in cgmodel.sidechain_positions:
                  parent_bead = bead_index - cgmodel.sidechain_length - 1
                 else:
                  parent_bead = bead_index - 1
              return(parent_bead)

def build_system(cgmodel):

        sigma = cgmodel.sigma.in_units_of(unit.nanometer)._value
        charge = cgmodel.charge._value
        epsilon = cgmodel.epsilon.in_units_of(unit.kilojoule_per_mole)._value
        bond_length = cgmodel.bond_length.in_units_of(unit.nanometer)._value

        # Create system
        system = mm.System()
        nonbonded_force = mm.NonbondedForce()
        bead_index = 0

        # Create nonbonded forces
        for monomer in range(cgmodel.polymer_length):
          for backbone_bead in range(cgmodel.backbone_length):
            system.addParticle(cgmodel.mass)
            nonbonded_force.addParticle(charge,sigma,epsilon)
            if backbone_bead in cgmodel.sidechain_positions:
              for sidechain_bead in range(cgmodel.sidechain_length):
                system.addParticle(cgmodel.mass)
                nonbonded_force.addParticle(charge,sigma,epsilon)

        # Create harmonic (bond potential) forces
        bond_list = cgmodel.get_bond_list()
        new_bond_list = []
        bead_index = 1
        for bond in bond_list:
              new_bond = [-1,-1]
              new_bond[0],new_bond[1] = bond[0]-1,bond[1]-1
              new_bond_list.append(new_bond)
              force = mm.HarmonicBondForce()
              if not cgmodel.constrain_bonds:
               force.addBond(new_bond[0],new_bond[1],bond_length,cgmodel.bond_force_constant)
              if cgmodel.constrain_bonds:
               system.addConstraint(new_bond[0],new_bond[1], bond_length)
               force.addBond(new_bond[0],new_bond[1],bond_length,0.0)
              system.addForce(force)

        
        for interaction in cgmodel.get_nonbonded_interaction_list():
               nonbonded_force.addException(interaction[0],interaction[1],cgmodel.charge,cgmodel.sigma,0.0)

        system.addForce(nonbonded_force)
        return(system)

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

        bond_force_constant: Bond force constant for all beads that are bonded,
        ( float )
        default = 9.9e5 kJ/mol/A^2


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
        bond_force_constant
        bb_bond_length
        bs_bond_length
        ss_bond_length
        charge
        num_beads
        positions
        system

        Notes
        -----
        
        """

        # Built in class attributes
        _BUILT_IN_REGIONS = ('polymer_length','backbone_length','sidechain_length','sidechain_positions','mass','sigma','epsilon','bond_length','bond_force_constant','bs_bond_length','bb_bond_length','ss_bond_length','charge','num_beads','positions','system','topology','constrain_bonds','bond_list','nonbonded_interactions')

        def __init__(self, positions = None, polymer_length = 12, backbone_length = 1, sidechain_length = 1, sidechain_positions = [0], mass = 12.0 * unit.amu, sigma = 8.4 * unit.angstrom, epsilon = 0.5 * unit.kilocalorie_per_mole, bond_length = 1.0 * unit.angstrom, bond_force_constant = 9.9e5, bb_bond_length = 1.0 * unit.angstrom, bs_bond_length = 1.0 * unit.angstrom, ss_bond_length = 1.0 * unit.angstrom, charge = 0.0 * unit.elementary_charge,constrain_bonds = False):

          """
          Initialize variables that were passed as input
          """

          self.polymer_length = polymer_length
          self.backbone_length = backbone_length
          self.sidechain_length = sidechain_length
          self.sidechain_positions = sidechain_positions
          self.num_beads = polymer_length * ( backbone_length + sidechain_length )
          self.mass = mass
          self.sigma = sigma
          self.epsilon = epsilon
          self.bond_length = bond_length
          self.bond_force_constant = bond_force_constant
          self.bb_bond_length = bb_bond_length
          self.bs_bond_length = bs_bond_length
          self.ss_bond_length = ss_bond_length
          self.charge = charge
          self.bond_list = self.get_bond_list()
          self.nonbonded_interactions = self.get_nonbonded_interaction_list()
          self.constrain_bonds = constrain_bonds

          """
          Initialize new (coarse grained) particle types:
          """

          """
          Make a list of coarse grained particle masses:
          """
          list_of_masses = get_particle_masses(self)

          add_new_elements(self,list_of_masses)

          self.system = build_system(self)

          self.positions = util.random_positions(self) 

          """
          Initialize attributes of our coarse grained model.
          """

        def get_bond_list(self):
          bond_list = []
          bead_index = 1
          for monomer in range(self.polymer_length):
            for backbone_bead in range(1,self.backbone_length+1):

             parent_index = get_parent_bead(self,bead_index,backbone_bead,sidechain_bead=False)
             if parent_index < -1:
              print("Error: identifying parent index incorrectly when assigning a bond.")
              print("The bead index is: "+str(bead_index))
              print("The parent index is: "+str(parent_index))
              print("The backbone index is: "+str(backbone_bead))
              exit() 
             if parent_index != -1:
              if parent_index < bead_index:
               bond_list.append([parent_index,bead_index])
              else:
               bond_list.append([bead_index,parent_index])
             bead_index = bead_index + 1
             
             if backbone_bead > 1:
              if backbone_bead-1 in self.sidechain_positions:
               for sidechain_bead in range(self.sidechain_length):
                 parent_index = get_parent_bead(self,bead_index,backbone_bead,sidechain_bead=True)
                 if parent_index < bead_index:
                  bond_list.append([parent_index,bead_index])
                 else:
                  bond_list.append([bead_index,parent_index])
                 bead_index = bead_index + 1

          return(bond_list)

        def get_nonbonded_interaction_list(cgmodel):
             interaction_list = []
             bond_list = [[bond[0]-1,bond[1]-1] for bond in cgmodel.get_bond_list()]
             for particle_1 in range(cgmodel.num_beads):
               for particle_2 in range(cgmodel.num_beads):
                 if particle_1 != particle_2 and abs(particle_1 - particle_2) >= 3:
                   if [particle_1,particle_2] not in bond_list and [particle_2,particle_1] not in bond_list:
                     if [particle_1,particle_2] not in interaction_list:
                       if [particle_2,particle_1] not in interaction_list:
                         interaction_list.append([particle_1,particle_2])
                     if [particle_2,particle_1] not in interaction_list:
                       if [particle_1,particle_2] not in interaction_list:
                         interaction_list.append([particle_2,particle_1])
             return(interaction_list)


        def get_dihedral_angles(self):
          bead_index = 0
          backbone_bead_indices = []
          dihedrals = []
          for monomer in range(self.polymer_length):
           for backbone_bead in range(self.backbone_length):
            backbone_bead_indices.append(bead_index)
            if bead_index != 0:
              bead_index = bead_index + 1
            if backbone_bead in self.sidechain_positions:
             for sidechain_bead in range(self.sidechain_length):
               bead_index = bead_index + 1

          for index in range(4,len(backbone_bead_indices)):

            dihedrals.append(np.array(backbone_bead_indices[index-4:index]))

          dihedrals = np.array([dihedral for dihedral in dihedrals])
          return(dihedrals)
