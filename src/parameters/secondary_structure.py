import os, subprocess
import numpy as np
from simtk import unit
from statistics import mean
from scipy.stats import linregress
from scipy import spatial
#from spatial.transform import Rotation as R
from foldamers.src.utilities.util import *
from foldamers.src.utilities.iotools import write_pdbfile_without_topology

def fraction_native_contacts(cgmodel,positions,native_structure,cutoff_distance=None):
        """
        """
        cutoff_distance = 1.1 * cgmodel.get_sigma(0)

        nonbonded_interaction_list = cgmodel.nonbonded_interaction_list
        #print("There are "+str(len(nonbonded_interaction_list))+" total nonbonded interaction possibilities.")
        native_structure_distances = distances(nonbonded_interaction_list,native_structure)
        current_structure_distances = distances(nonbonded_interaction_list,positions)
        native_distances = []
        native_interaction_list = []
        
        for interaction in range(len(nonbonded_interaction_list)):
          if native_structure_distances[interaction].__lt__(cutoff_distance):
            native_distances.append(native_structure_distances[interaction])
            native_interaction_list.append(interaction)
        total_native_interactions = len(native_interaction_list)

        current_distances = []
        current_structure_interaction_list = []
        for interaction in range(len(nonbonded_interaction_list)):
          if interaction in native_interaction_list:
            if current_structure_distances[interaction].__lt__(cutoff_distance):
              current_distances.append(current_structure_distances[interaction])
              current_structure_interaction_list.append(interaction)
        current_structure_native_interactions = len(current_structure_interaction_list) 
        Q = current_structure_native_interactions / total_native_interactions
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

def get_helical_data(cgmodel):
        """
        """
        file = open("before_rotation.pdb","w")
        PDBFile.writeFile(cgmodel.topology,cgmodel.positions,file=file)
        positions = np.array([[float(i.in_units_of(unit.angstrom)._value) for i in position] for position in cgmodel.positions])
        # 1) Get the backbone particle positions
        backbone_positions = []
        for particle in range(len(cgmodel.positions)):
          if cgmodel.get_particle_type(particle) == "backbone":
            backbone_positions.append(cgmodel.positions[particle])
        backbone_positions = np.array([[float(i.in_units_of(unit.angstrom)._value) for i in coord] for coord in backbone_positions])
        # 2) Project the backbone positions onto the (x,y) plane
        print(backbone_positions)
        xy_projected_positions = backbone_positions
        x_data = []
        y_data = []
        for position in xy_projected_positions:
          position[2] = 0.0
          x_data.append(position[0])
          y_data.append(position[1])
        # 3) Calculate the best fit line for these projected positions
        slope,intercept,r,p,std_err=linregress(np.array([x for x in x_data]),np.array([y for y in y_data]))
        # 4) Rotate the coordinate system so that this line is oriented along the x-axis
        # Calculate angle from z-axis:
        z_axis_angle = np.arctan(slope)
        z_axis_rotation_matrix = spatial.transform.Rotation.from_euler('z', z_axis_angle, degrees=False)
        x_oriented_positions = z_axis_rotation_matrix.apply(positions)
        # 5) Project the positions onto the (x,z) plane
        xz_projected_positions = backbone_positions
        x_data = []
        z_data = []
        for position_index in range(len(xz_projected_positions)):
          xz_projected_positions[position_index][1] = 0.0
          x_data.append(xz_projected_positions[position_index][0])
          z_data.append(xz_projected_positions[position_index][2])

        # 6) Calculate the best fit line for these projected positions
        slope,intercept,r,p,std_err=linregress(np.array([x for x in x_data]),np.array([z for z in z_data]))
        # 7) Rotate the coordinate system so that this line is oriented along the x-axis
        # Calculate angle from y-axis:
        y_axis_angle = np.arctan(slope)
        y_axis_rotation_matrix = spatial.transform.Rotation.from_euler('y', y_axis_angle, degrees=False)
        final_positions = y_axis_rotation_matrix.apply(x_oriented_positions)

        cgmodel.positions = unit.Quantity(final_positions,cgmodel.positions.unit)

        file = open("after_rotation.pdb","w")
        PDBFile.writeFile(cgmodel.topology,cgmodel.positions,file=file)

        # 8) Using these transformed coordinates, calculate the helical parameters for this structure:

        # radius
        axis_distances = []
        rotations = 0.0
        for position in final_positions:
          axis_distance = distance(unit.Quantity([float(position[0]),0.0,0.0],cgmodel.positions.unit),unit.Quantity(position,cgmodel.positions.unit))
          axis_distances.append(axis_distance)
          if len(axis_distances) > 1:
            rotation = np.arctan(position[1]/position[2]) - last_angle
            last_angle = rotation
            rotations = rotations + rotation
          else:
            rotation = np.arctan(position[1]/position[2])
            last_angle = rotation
            rotations = rotations + rotation

            
        print(axis_distances)
        exit()
        radius = mean(np.array([float(dist.in_units_of(unit.angstrom)._value) for dist in axis_distances]))
        particles_per_turn = float(cgmodel.polymer_length/(rotations/6.28))

        # pitch
        #
        # Shift all coordinates so that the first backbone atom has x=0

        shift = - final_positions[0][0]

        axis_deltas = []
        for position in final_positions:
          position[0] = position[0] + shift
          if abs(position[0] - final_positions[0][0]) > 0:
            axis_deltas.append(float(position[0]-final_positions[0][0]))
        average_delta = mean(axis_deltas)
        pitch = average_delta * particles_per_turn

        return(radius,pitch,particles_per_turn)
        

