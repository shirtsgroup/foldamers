import os
from simtk import unit
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
from foldamers.cg_model.cgmodel import CGModel
from foldamers.ensembles.cluster import *
from cg_openmm.simulation.rep_exch import *

# Coarse grained model settings
polymer_length = 24
backbone_lengths = [1]
sidechain_lengths = [1]
sidechain_positions = [0]
include_bond_forces = True
include_bond_angle_forces = True
include_nonbonded_forces = True
include_torsion_forces = False
constrain_bonds = False
random_positions = True

# Replica exchange simulation settings
total_simulation_time = 0.5 * unit.nanosecond
simulation_time_step = 10.0 * unit.femtosecond
total_steps = int(np.floor(total_simulation_time / simulation_time_step))
number_replicas = 48
min_temp = 50.0 * unit.kelvin
max_temp = 600.0 * unit.kelvin
exchange_frequency = 100  # Number of steps between exchange attempts

# Bond definitions
bond_length = 7.5 * unit.angstrom
bond_lengths = {
    "bb_bb_bond_length": bond_length,
    "bb_sc_bond_length": bond_length,
    "sc_sc_bond_length": bond_length,
}
bond_force_constant = 0 * unit.kilocalorie_per_mole / unit.nanometer / unit.nanometer
bond_force_constants = {
    "bb_bb_bond_k": bond_force_constant,
    "bb_sc_bond_k": bond_force_constant,
    "sc_sc_bond_k": bond_force_constant,
}

# Particle definitions
mass = 100.0 * unit.amu
masses = {"backbone_bead_masses": mass, "sidechain_bead_masses": mass}
r_min = 1.5 * bond_length  # Lennard-Jones potential r_min
# Factor of /(2.0**(1/6)) is applied to convert r_min to sigma
sigma = r_min / (2.0 ** (1.0 / 6.0))
sigmas = {"bb_sigma": sigma, "sc_sigma": sigma}
epsilon = 0.5 * unit.kilojoule_per_mole
epsilons = {"bb_eps": epsilon, "sc_eps": epsilon}

# Bond angle definitions
bond_angle_force_constant = 100 * unit.kilojoule_per_mole / unit.radian / unit.radian
bond_angle_force_constants = {
    "bb_bb_bb_angle_k": bond_angle_force_constant,
    "bb_bb_sc_angle_k": bond_angle_force_constant,
}
# OpenMM requires angle definitions in units of radians
bb_bb_bb_equil_bond_angle = 120.0 * unit.degrees
bb_bb_sc_equil_bond_angle = 120.0 * unit.degrees
equil_bond_angles = {
    "bb_bb_bb_angle_0": bb_bb_bb_equil_bond_angle,
    "bb_bb_sc_angle_0": bb_bb_sc_equil_bond_angle,
}

# Torsion angle definitions
torsion_force_constant = 20.0 * unit.kilojoule_per_mole
torsion_force_constants = {"bb_bb_bb_bb_torsion_k": torsion_force_constant}
# OpenMM requires angle definitions in units of radians
bb_bb_bb_bb_equil_torsion_angle = 78.0 * unit.degrees
bb_bb_bb_sc_equil_torsion_angle = 78.0 * unit.degrees
equil_torsion_angles = {"bb_bb_bb_bb_torsion_0": bb_bb_bb_bb_equil_torsion_angle}
torsion_periodicities = {"bb_bb_bb_bb_period": 3}

# Get initial positions from local file
positions = PDBFile("24mer_1b1s_initial_structure.pdb").getPositions()

# Build a coarse grained model
cgmodel = CGModel(
    polymer_length=polymer_length,
    backbone_lengths=backbone_lengths,
    sidechain_lengths=sidechain_lengths,
    sidechain_positions=sidechain_positions,
    masses=masses,
    sigmas=sigmas,
    epsilons=epsilons,
    bond_lengths=bond_lengths,
    bond_force_constants=bond_force_constants,
    bond_angle_force_constants=bond_angle_force_constants,
    torsion_force_constants=torsion_force_constants,
    equil_bond_angles=equil_bond_angles,
    equil_torsion_angles=equil_torsion_angles,
    torsion_periodicities=torsion_periodicities,
    include_nonbonded_forces=include_nonbonded_forces,
    include_bond_forces=include_bond_forces,
    include_bond_angle_forces=include_bond_angle_forces,
    include_torsion_forces=include_torsion_forces,
    constrain_bonds=constrain_bonds,
    positions=positions,
)


# Create list of trajectory files for clustering analysis	
pdb_file_list = []
for i in range(number_replicas):
    pdb_file_list.append("output/replica_%s.pdb" %(i+1))

# Set clustering parameters
n_clusters=2
frame_start=0
frame_stride=10
frame_end=-1

# Run KMeans clustering
centroid_positions, medoid_positions = get_cluster_centroid_positions(
    pdb_file_list=pdb_file_list,
    cgmodel=cgmodel,
    n_clusters=n_clusters,
    frame_start=frame_start,
    frame_stride=frame_stride,
    frame_end=-1)
