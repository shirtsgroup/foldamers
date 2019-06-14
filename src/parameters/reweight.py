#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as pyplot
# OpenMM utilities
import mdtraj as md
from simtk import unit
# foldamers utilities
from foldamers.src.cg_model.cgmodel import basic_cgmodel
from foldamers.src.ensembles.ens_build import 

def get_expectation_value(data_set):
        """
        Evaluate the dimensionless 

        Parameters
        ----------

        cgmodel: CGModel() class object.

        Returns
        -------

        cgmodel: CGModel() class object.

        """

        # Set variable model settings
        base_sigma = cgmodel.sigmas['bb_bb_sigma'] # Lennard-Jones interaction distance
        base_epsilon = cgmodel.epsilons['bb_bb_epsilon'] # Lennard-Jones interaction strength
        sigma_list = [(base_sigma).__add__(i * base_sigma.unit) for i in [ j * 0.2 for j in range(-2,3,1)]]
        epsilon_list = [(base_epsilon).__add__(i * base_epsilon.unit) for i in [ j * 0.2 for j in range(-1,3,1)]]
        sigma_epsilon_list = np.zeros((len(sigma_list),len(epsilon_list)))

        for sigma_index in range(len(sigma_list)):
          for epsilon_index in range(len(epsilon_list)):
            sigma = sigma_list[sigma_index]
            epsilon = epsilon_list[epsilon_index]
            print("Evaluating the energy for a model with:")
            print("sigma="+str(sigma)+" and epsilon="+str(epsilon))
            # Build a coarse grained model
            cgmodel = basic_cgmodel(polymer_length=polymer_length, backbone_length=backbone_length, sidechain_length=sidechain_length, sidechain_positions=sidechain_positions, mass=mass, sigma=sigma, epsilon=epsilon, bond_length=bond_length)
            # Get the average energy for an unfolded (high-energy) ensemble of structures built with these model settings
            
            energy = min(energy['potential_energy'])
            sigma_epsilon_list[sigma_index][epsilon_index] = minimum_energy


        return(cgmodel)
