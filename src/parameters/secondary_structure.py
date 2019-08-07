import numpy as np
from foldamers.src.utilities.util import distances

def fraction_native_contacts(cgmodel,native_structure,cutoff_distance=None):
        """
        """
        cutoff_distance = 1.2 * cgmodel.get_sigma(0)

        nonbonded_interaction_list = cgmodel.nonbonded_interaction_list
        #print("There are "+str(len(nonbonded_interaction_list))+" total nonbonded interaction possibilities.")
        native_structure_distances = distances(nonbonded_interaction_list,native_structure)
        current_structure_distances = distances(nonbonded_interaction_list,cgmodel.positions)
        native_distances = []
        native_interaction_list = []
        
        for interaction in range(len(nonbonded_interaction_list)):
          if native_structure_distances[interaction].__lt__(cutoff_distance):
            #print(native_structure_distances[interaction])
            #print(cutoff_distance.in_units_of(native_structure_distances[interaction].unit))
            native_distances.append(native_structure_distances[interaction])
            native_interaction_list.append(interaction)
        total_native_interactions = len(native_interaction_list)
        #print("There are "+str(total_native_interactions)+" total nonbonded interactions within the cutoff")
        #print("distance for the native structure.")
        #print("Their distances are: "+str([dist._value for dist in native_distances]))

        current_distances = []
        current_structure_interaction_list = []
        for interaction in range(len(nonbonded_interaction_list)):
          if interaction in native_interaction_list:
            if current_structure_distances[interaction].__lt__(cutoff_distance):
              current_distances.append(current_structure_distances[interaction])
              current_structure_interaction_list.append(interaction)
        current_structure_native_interactions = len(current_structure_interaction_list) 
        #print("There are "+str(current_structure_native_interactions)+" total 'native' interactions in the current structure.")
        #print("Their distances are: "+str([dist._value for dist in current_distances]))
        Q = current_structure_native_interactions / total_native_interactions
        #print(Q)
        return(Q)

#def 
