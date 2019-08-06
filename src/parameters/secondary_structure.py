import numpy as np
from foldamers.src.utilities.util import distances

def fraction_native_contacts(cgmodel,native_structure,cutoff_distance=None):
        """
        """
        if cutoff_distance == None:
          cutoff_distance = 1.5 * cgmodel.get_sigma(0)

        nonbonded_interaction_list = cgmodel.nonbonded_interaction_list
        native_structure_distances = distances(nonbonded_interaction_list,native_structure)
        current_structure_distances = distances(nonbonded_interaction_list,cgmodel.positions)
        native_interactions = 0
        for interaction in range(len(nonbonded_interaction_list)):
          if native_structure_distances[interaction] < cutoff_distance:
            native_interactions = native_interactions + 1
        total_native_interactions = 0
        for interaction in range(len(nonbonded_interaction_list)):
          if native_structure_distances[interaction] < cutoff_distance:
            if current_structure_distances[interaction] < cutoff_distance:
              total_native_interactions = total_native_interactions + 1
        Q = total_native_interactions / native_interactions
        return(Q)

#def 
