## No default python environment

import os, sys, timeit
from io import StringIO
import numpy as np
import math, random
import matplotlib.pyplot as pyplot
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
total_simulation_time = 100.0 * unit.picosecond # Units = picoseconds

# Model settings
backbone_length = 1 # Number of backbone beads
sidechain_length = 1 # Number of sidechain beads
sidechain_positions = [0] # Index of backbone bead on which the side chains are placed
polymer_length = 12 # Number of monomers in the polymer
mass = 12.0 * unit.amu # Mass of beads
sigma = 8.4 * unit.angstrom # Lennard-Jones interaction distance
bond_length = 1.0 * unit.angstrom # bond length
bond_force_constant = 9.9e5 # Units = kJ/mol/A^2
constrain_bonds = False
epsilon = 0.5 * unit.kilocalorie_per_mole # Lennard-Jones interaction strength
charge = 0.0 * unit.elementary_charge # Charge of beads

sigma = sigma.in_units_of(unit.nanometer)._value
charge = charge._value
epsilon = epsilon.in_units_of(unit.kilojoule_per_mole)._value
bond_length = bond_length.in_units_of(unit.nanometer)._value

def add_new_elements():
        elem.Element(117,'cg-backbone','CG1',mass)
        elem.Element(118,'cg-sidechain','CG2',mass)
        return

def build_system():
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
cgmodel = CGModel()
system = build_system()
pdb_file = "test.pdb"
write_pdbfile(cgmodel,pdb_file)
pdb_mm_obj = PDBFile(pdb_file)
topology = pdb_mm_obj.getTopology()
simulation = build_mm_simulation(topology,system,cgmodel.positions,temperature=temperature,simulation_time_step=simulation_time_step,total_simulation_time=simulation_time_step*print_frequency,output_pdb="dihedrals.pdb",output_data="dihedrals_sim_data.dat",print_frequency=print_frequency)
energies=[]
number_stages = round((total_simulation_time._value/simulation_time_step._value)/print_frequency)
print("Number of simulation stages is: "+str(number_stages))
for stage in range(number_stages):
 simulation.step(print_frequency)
 energy = round(simulation.context.getState(getEnergy=True).getPotentialEnergy()._value,2)
 energies.append(energy)
trajectory = md.load("dihedrals.pdb")
dihedral_list = get_dihedral_angles()
dihedrals = []
for dihedral in dihedral_list:
 angles = md.compute_dihedrals(trajectory,[dihedral])
 dihedrals.append(angles)
figure = pyplot.figure(0)
colors = pyplot.cm.rainbow(np.linspace(0, 1, len(dihedral_list)))
for dihedral, c in zip(dihedrals, colors):
  pyplot.scatter(dihedral, energies, color=c, figure=figure)
pyplot.xlabel("Dihedral Angle (Degrees)")
pyplot.ylabel("Potential Energy (kJ/mol)")
pyplot.title("Dihedral distribution data for simulation of 1,1-CG model")
pyplot.savefig("/CG-11-distribution.png")
pyplot.legend(dihedral_list)
pyplot.show()
pyplot.close()
exit()
