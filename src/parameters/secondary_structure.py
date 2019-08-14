import os, subprocess
import numpy as np
from foldamers.src.utilities.util import distances
from foldamers.src.utilities.iotools import write_pdbfile_without_topology

def fraction_native_contacts(cgmodel,positions,native_structure,cutoff_distance=None):
        """
        """
        cutoff_distance = 1.2 * cgmodel.get_sigma(0)

        nonbonded_interaction_list = cgmodel.nonbonded_interaction_list
        #print("There are "+str(len(nonbonded_interaction_list))+" total nonbonded interaction possibilities.")
        native_structure_distances = distances(nonbonded_interaction_list,native_structure)
        current_structure_distances = distances(nonbonded_interaction_list,positions)
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

def get_helical_parameters(cgmodel):
        """
        """
        helios_path = str(str(os.path.realpath(__file__).split('/secondary_structure')[0])+str("/helios.o"))
        write_pdbfile_without_topology(cgmodel,"temp_pitch.pdb")
        kHelix_run_file = "run_kHelix.sh"
        file = open(kHelix_run_file,"w")
        file.write('#!/bin/bash\n')
        file.write('\n')
        file.write('cat > input << EOF\n')
        file.write('inputhelix $1\n')
        file.write('helixout_name kHelix.out\n')
        file.write('coord_type 1\n')
        file.write('num_grid 360\n')
        file.write('natoms '+str(round(cgmodel.num_beads/2))+'\n')
        file.write('nframes 1\n')
        file.write('grid_phi_beg 0\n')
        file.write('grid_phi_end 180\n')
        file.write('grid_theta_beg 0\n')
        file.write('grid_theta_end 180\n')
        file.write('helix_atom_names X1\n')
        file.write('print_to_plot 1\n')
        file.write('EOF\n')
        file.write(str(helios_path)+' input\n')
        #file.write('done\n')        
        file.close()
        subprocess.run([str(str(os.getcwd())+"/"+str(kHelix_run_file)),"temp_pitch.pdb",">","helios_output"])
        #os.remove("helios_output")
        file = open("kHelix.out",mode="r")
        output = file.readlines()
        line_index = 1
        for line in output:
          if line_index == 43:
            radius = line.split()[3]
            pitch = line.split()[4]
            sweep = float(line.split()[5])
            monomers_per_turn = cgmodel.polymer_length/(sweep/360.0)
            break
          line_index = line_index + 1
          
        return(pitch,radius,monomers_per_turn)

