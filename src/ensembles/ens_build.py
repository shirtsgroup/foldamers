#!/usr/bin/python

import os
from simtk import unit
from foldamers.src.cg_model.cgmodel import basic_cgmodel
from cg_openmm.src.build.cg_build import get_mm_energy
from foldamers.src.utilities.util import random_positions



def get_ensemble(cgmodel,ensemble_size=100,high_energy=False,low_energy=False):
        """
        Given a coarse grained model, this function generates an ensemble of high energy configurations and, by default, saves this ensemble to the foldamers/ensembles database for future reference/use, if a high-energy ensemble with these settings does not already exist.

        Parameters
        ----------

        :param cgmodel: CGModel() class object.
        :type cgmodel: class

        :param ensemble_size: Number of structures to generate for this ensemble, default = 100
        :type ensemble_size: integer

        :param high_energy: If set to 'True', this function will generate an ensemble of high-energy structures, default = False
        :type high_energy: Logical

        :param low_energy: If set to 'True', this function will generate an ensemble of low-energy structures, default = False
        :type low_energy: Logical

        Returns
        -------

        ensemble: List( positions( np.array( float * simtk.unit ( shape = num_beads x 3 ) ) )
                  A list of the positions for all members in the high_energy ensemble.

        """
        if high_energy and low_energy:
          print("ERROR: Both 'high_energy' and 'low_energy' ensembles were requested in 'get_ensemble()'.  Please set only one of these variables to 'True', and call the function again.")
          exit()
        if low_energy:
          print("Generating an ensemble of "+str(ensemble_size)+" low energy configurations.")
        if high_energy:
          print("Generating an ensemble of "+str(ensemble_size)+" high energy configurations.")
        if not high_energy and not low_energy:
          print("Generating an ensemble of "+str(ensemble_size)+" configurations.")

        ensemble = []
        for member in range(ensemble_size):

          if high_energy:
            positions = random_positions(cgmodel,high_energy=True)
            
          if low_energy:
            positions = random_positions(cgmodel,low_energy=True)

          if not high_energy and not low_energy:
            positions = random_positions(cgmodel)

          ensemble.append(positions)
        
        return(ensemble)

def z_score(topology,system,nonnative_ensemble,native_structure):
        """
        Given an ensemble of nonnative structures, and a low-energy ("native") structure, this subroutine will calculate the Z-score.

        Parameters
        ----------

        nonnative_ensemble: List( positions( np.array( float * simtk.unit ( shape = num_beads x 3 ) ) )
                  A list of the positions for all members in the high_energy ensemble.

        native_structure: positions( np.array( float * simtk.unit ( shape = num_beads x 3 ) )
                          The positions for a low energy structure.

        """

        nonnative_ensemble_energies = []

        for pose in nonnative_ensemble:
         
          energy = get_mm_energy(topology,system,pose)
          nonnative_ensemble_energies.append(energy)

        average_nonnative_energy = statistics.mean(nonnative_ensemble_energies)

        stdev_nonnative_energy = statistics.stdev(nonnative_ensemble_energies)

        native_energy = get_mm_energy(topology,system,native_structure)

        z_score = ( average_nonnative_energy - native_energy ) / stdev_nonnative_energy

        return(z_score)
