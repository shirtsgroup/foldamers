## No default python environment

import os, sys, timeit
from io import StringIO
import numpy as np
import math, random
import matplotlib.pyplot as pyplot
import statistics
from statistics import mean
# OpenMM utilities
import simtk.openmm.app.element as elem
from simtk.openmm.app.pdbfile import PDBFile
import mdtraj as md
from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk import unit
from simtk.openmm.vec3 import Vec3
# foldamers utilities
from foldamers.src.cg_model.cgmodel import CGModel
from foldamers.src.utilities.iotools import *
from cg_openmm.src.cg_mm_tools.cg_openmm import *

box_size = 10.00 * unit.nanometer # box width
cutoff = box_size / 2.0 * 0.99
simulation_time_step = 0.002 * unit.picosecond # Units = picoseconds
temperature = 600.0 * unit.kelvin
print_frequency = 10 # Number of steps to skip when printing output
total_simulation_time = 1.0 * unit.picosecond # Units = picoseconds

# Model settings
backbone_length = 1 # Number of backbone beads
sidechain_length = 1 # Number of sidechain beads
sidechain_positions = [0] # Index of backbone bead on which the side chains are placed
polymer_length = 12 # Number of monomers in the polymer
mass = 12.0 * unit.amu # Mass of beads
bond_length = 1.0 * unit.angstrom # bond length
bond_force_constant = 9.9e5 # Units = kJ/mol/A^2
constrain_bonds = False
charge = 0.0 * unit.elementary_charge # Charge of beads

ensemble_size = 100
number_sigma = 20
number_epsilon = 20
min_sigma = 5.0
min_epsilon = 0.1
sigma_step = 0.5
epsilon_step = 0.2

charge = charge._value
bond_length = bond_length.in_units_of(unit.nanometer)._value

def add_new_elements():
        elem.Element(117,'cg-backbone','CG1',mass)
        elem.Element(118,'cg-sidechain','CG2',mass)
        return

def build_system(charge,sigma,epsilon):
        # Create system
        system = mm.System()
        nonbonded_force = mm.NonbondedForce()
        bead_index = 0
        for monomer in range(polymer_length):
          for backbone_bead in range(backbone_length):
            system.addParticle(mass)
            nonbonded_force.addParticle(charge,sigma,epsilon)
            if monomer != 0:
             bead_index = bead_index + 1

             if backbone_bead == 0:
              force = mm.HarmonicBondForce()
              force.addBond(bead_index-sidechain_length-1, bead_index, bond_length,bond_force_constant)
              system.addForce(force)
              nonbonded_force.addException(bead_index-sidechain_length-1,bead_index,charge,bond_length,epsilon=0.0)

              if constrain_bonds:
               system.addConstraint(bead_index-sidechain_length-1, bead_index, bond_length)

             if backbone_bead != 0:
              force = mm.HarmonicBondForce()
              force.addBond(bead_index-1, bead_index, bond_length,bond_force_constant)
              system.addForce(force)
              nonbonded_force.addException(bead_index-1, bead_index,charge,bond_length,epsilon=0.0)

              if constrain_bonds:
               system.addConstraint(bead_index-1, bead_index, bond_length)

            if backbone_bead in sidechain_positions:
              for sidechain in range(sidechain_length):
                system.addParticle(mass)
                nonbonded_force.addParticle(charge,sigma,epsilon)
                bead_index = bead_index + 1

                force = mm.HarmonicBondForce()
                force.addBond(bead_index-1, bead_index, bond_length,bond_force_constant)
                system.addForce(force)
                nonbonded_force.addException(bead_index-1, bead_index,charge,bond_length,epsilon=0.0)

                if constrain_bonds:
                  system.addConstraint(bead_index,bead_index-1,bond_length)
  
        system.addForce(nonbonded_force)

        return(system)

def get_dihedral_angles():
  bead_index = 0
  backbone_bead_indices = []
  dihedrals = []
  for monomer in range(polymer_length):
    for backbone_bead in range(backbone_length):
      backbone_bead_indices.append(bead_index)
      if bead_index != 0:
        bead_index = bead_index + 1
      for sidechain_bead in range(sidechain_length):
        bead_index = bead_index + 1

  for index in range(4,len(backbone_bead_indices)):
    
    dihedrals.append(np.array(backbone_bead_indices[index-4:index]))

  dihedrals = np.array([dihedral for dihedral in dihedrals])
  return(dihedrals)

add_new_elements()

sigma = [ float(min_sigma + i * sigma_step) * unit.angstrom for i in range(number_sigma)] # Lennard-Jones interaction distance
epsilon = [ float(min_epsilon + i * epsilon_step) * unit.kilojoule_per_mole for i in range(number_epsilon)] # Lennard-Jones interaction strength
potential_energies = np.zeros([len(sigma),len(epsilon)])

for sig_index in range(len(sigma)):

 for eps_index in range(len(epsilon)):

   ensemble_energies=[]
   for pose in range(ensemble_size):

     cgmodel = CGModel()
     pdb_file = "test.pdb"
     write_pdbfile(cgmodel,pdb_file)
     pdb_mm_obj = PDBFile(pdb_file)
     topology = pdb_mm_obj.getTopology()
     sig = sigma[sig_index].in_units_of(unit.nanometer)._value
     eps = epsilon[eps_index].in_units_of(unit.kilojoule_per_mole)._value
     system = build_system(charge,sig,eps)
     simulation = build_mm_simulation(topology,system,cgmodel.positions,temperature=temperature,simulation_time_step=simulation_time_step,total_simulation_time=simulation_time_step,output_pdb="opt_sig_eps.pdb",output_data="opt_sig_eps.dat",print_frequency=print_frequency)

     energy = round(simulation.context.getState(getEnergy=True).getPotentialEnergy()._value,2)
     ensemble_energies.append(energy)

   potential_energies[sig_index][eps_index] = mean(ensemble_energies)

sigma = [sigma[i]._value for i in range(len(sigma))]
epsilon = [epsilon[i]._value for i in range(len(epsilon))]
x=np.unique(sigma)
y=np.unique(epsilon)
z=np.array(potential_energies)
X,Y = np.meshgrid(x,y)

Z=z.reshape(len(y),len(x))

figure_index = 1
figure = pyplot.figure(figure_index)
pyplot.ylabel('Epsilon ( kJ/mol )')
pyplot.xlabel('Sigma ( Angstroms )')
heatmap = pyplot.pcolor(X,Y,Z)
cbar = pyplot.colorbar(heatmap)
cbar.set_label('Potential Energy ( kJ/mol )', rotation=270)
pyplot.savefig(str("sigma_epsilon_heatmap.png"))
pyplot.show()
pyplot.close()


exit()
