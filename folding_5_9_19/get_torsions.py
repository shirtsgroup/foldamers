#!/usr/bin/python

# Import external Python packages
import numpy as np
import matplotlib.pyplot as pyplot
from simtk.openmm.app.pdbfile import PDBFile
from simtk import unit
import mdtraj as md

# foldamers utilities
from foldamers.src.cg_model.cgmodel import CGModel
from foldamers.src.utilities.util import *
from foldamers.src.utilities.iotools import *
from cg_openmm.src.cg_mm_tools.cg_openmm import *

# Simulation settings
input_coordinates = "test.pdb"
output_coordinates = "simulation.pdb"
output_data = "simulation.dat"
box_size = 10.00 * unit.nanometer # box width
cutoff = box_size / 2.0 * 0.99
simulation_time_step = 0.5 * unit.femtosecond # Units = picoseconds
temperature = 300.0 * unit.kelvin
print_frequency = 10 # Number of steps to skip when printing output
total_simulation_time = 1.0 * unit.picosecond # Units = picoseconds

# Model settings
backbone_length = 1 # Number of backbone beads
sidechain_length = 1 # Number of sidechain beads
sidechain_positions = [0] # Index of backbone bead(s) where sidechains are attached
polymer_length = 12 # Number of monomers in the polymer
mass = 12.0 * unit.amu # Mass of beads
sigma = 1.0 * unit.angstrom # Lennard-Jones interaction distance
bond_length = 1.0 * unit.angstrom # bond length
bond_force_constant = 0 # Units = kJ/mol/A^2
constrain_bonds = True
epsilon = 0.1 * unit.kilocalorie_per_mole # Lennard-Jones interaction strength
charge = 0.0 * unit.elementary_charge # Charge of beads

# Read coordinates and build a coarse-grained model
pdb_mm_obj = PDBFile(input_coordinates)
positions = pdb_mm_obj.getPositions()
cgmodel = CGModel(positions=positions,polymer_length=polymer_length,backbone_length=backbone_length, sidechain_length=sidechain_length, sidechain_positions = sidechain_positions, mass = mass, sigma = sigma, epsilon = epsilon, bond_length = bond_length, bond_force_constant = bond_force_constant, charge = charge,constrain_bonds=constrain_bonds)
topology = pdb_mm_obj.getTopology()

# Setup and run OpenMM simulation
simulation = build_mm_simulation(topology,cgmodel.system,cgmodel.positions,temperature=temperature,simulation_time_step=simulation_time_step,total_simulation_time=simulation_time_step*print_frequency,output_pdb=output_coordinates,output_data=output_data,print_frequency=print_frequency)
number_steps = round(total_simulation_time.__div__(simulation_time_step))
print("Running a simulation for: "+str(number_steps)+" steps.")
simulation.step(number_steps)
print("Simulation run successful.")

# Analyze simulation output
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
