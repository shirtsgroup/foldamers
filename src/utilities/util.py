
# =============================================================================================
# 1) PYTHON PACKAGE IMPORTS
# =============================================================================================

# System packages
import numpy as np
import math, random
from simtk import unit
from simtk.openmm.app.pdbfile import PDBFile
from cg_openmm.src.cg_mm_tools.cg_openmm import minimize_structure
from foldamers.src.cg_model.cgmodel import *
from foldamers.src.utilities.iotools import *

# =============================================================================================
# 2) ENVIRONMENT/JOB SETTINGS
# =============================================================================================

def random_sign(number):
        """
        Returns 'number' with a random sign.

        Parameters
        ----------

        number: float

        Returns
        -------

        number

        """

        # 0 = negative, 1 = positive
        random_int = random.randint(0,1)
        if random_int == 0: number = -1.0 * number

        return(number)

def first_bead(positions):
        """
        Determine if we have any particles in 'positions'
        Parameters
        ----------
        positions: Positions for all beads in the coarse-grained model.
        ( np.array( float * unit ( shape = num_beads x 3 ) ) )
        Returns
        -------
        first_bead: Logical variable stating if this is the first particle.
        """

        first_bead = True
    
        if str(positions.shape) == '(2,)':
          return(first_bead)
        else:
          for value in positions._value:
            if any(i != 0. for i in value):
              first_bead = False

        return(first_bead)


def get_move(trial_coordinates,move_direction,distance,bond_length,finish_bond=False):
        """
        Given a 'move_direction', a current distance, and a
        target 'bond_length' ( Index denoting x,y,z Cartesian 
        direction), update the coordinates for the particle.

        Parameters
        ----------

        trial_coordinates: positions for a particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        move_direction: Cartesian direction in which we will
        attempt a particle placement, where: x=0, y=1, z=2. 
        ( integer )

        distance: Current distance from parent particle
        ( float * simtk.unit.distance )

        bond_length: Target bond_length for particle placement.
        ( float * simtk.unit.distance )

        finish_bond: Logical variable determining how we will
        update the coordinates for this particle.

        Returns
        -------

        trial_coordinates: Updated positions for the particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        """

        if distance.__gt__(bond_length):
          print("ERROR: The particle distance is larger than the bond length.")
          exit()

        # Determine the 'max_step_size' as the square root of the difference
        # between 'bond_length' and 'distance'
#        print("Bond length is "+str(bond_length))
#        print("Distance is "+str(distance))
        max_step_size = bond_length.__pow__(2.0).__sub__(distance.__pow__(2.0)).sqrt()

        # Add a random sign to 'max_step_size', to randomize our particle placement.
        sign_index = random.randint(0,1)
        if sign_index == 0:
          max_step_size = max_step_size.__neg__()


        # If we are "finishing the bond", then the "step" size is
        # the length (in the direction 'move_direction') required
        # so that 'distance' == 'bond_length'
        if finish_bond:

          # Calculate the step size as 'step' = sqrt('difference')
          step = max_step_size._value

        # If we aren't "finishing the bond", then the "step" size
        # is a random float in the range 0.0
        if not finish_bond:

          step = random.uniform(0.0,max_step_size._value)

        # Add this 'step' to the existing coordinates
        # print("The trial coordinates are: "+str(trial_coordinates))
        # print("The step size is: "+str(step))
        trial_coordinates[move_direction] = trial_coordinates[move_direction].__add__(unit.Quantity(step,trial_coordinates.unit))
        # print("The trial coordinates are: "+str(trial_coordinates))

        return(trial_coordinates)

def attempt_lattice_move(parent_coordinates,bond_length,move_direction_list):
        """
        Given a set of cartesian coordinates, assign a new particle
        a distance of 'bond_length' away in a random direction.

        Parameters
        ----------

        parent_coordinates: Positions for a single particle,
        away from which we will place a new particle a distance
        of 'bond_length' away.
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        trial_coordinates: Positions for a new trial particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        """

        """ 'dist' tracks the distance between the new trial particle
        and the parent particle """

        units = bond_length.unit
        dist = unit.Quantity(0.0,units)

        """ 'move_direction_list' tracks the Cartesian
        directions in which we have attempted a particle placement,
        where: x=0,y=1,z=2. """

        # Assign the parent coordinates as the initial coordinates for a trial particle

        trial_coordinates = parent_coordinates.__deepcopy__(memo={})
        ref = parent_coordinates.__deepcopy__(memo={})

        move_direction = random.randint(1,3)
        move_sign = random_sign(move_direction)
        while move_sign in move_direction_list:
         move_direction = random.randint(1,3)
         move_sign = random_sign(move_direction)
        move_direction_list.append(move_sign)
        if move_sign < 0:
         move_sign = -bond_length
        else:
         move_sign = bond_length
#        print(move_direction)
#        print(trial_coordinates[move_direction-1])
#        print(move_sign)
        trial_coordinates[move_direction-1] = trial_coordinates[move_direction-1].__add__(move_sign)

        dist = distance(ref,trial_coordinates)

        if round(dist._value,4) < round(bond_length._value,4):

           print("Error: particles are being placed at a distance different from the bond length")
           print("Bond length is: "+str(bond_length))
           print("The particle distance is: "+str(dist))
           print(ref)
           print(trial_coordinates)
           exit()

        return(trial_coordinates,move_direction_list)

def attempt_move(parent_coordinates,bond_length):
        """
        Given a set of cartesian coordinates, assign a new particle
        a distance of 'bond_length' away in a random direction.

        Parameters
        ----------

        parent_coordinates: Positions for a single particle,
        away from which we will place a new particle a distance
        of 'bond_length' away.
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        trial_coordinates: Positions for a new trial particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        """

        """ 'dist' tracks the distance between the new trial particle
        and the parent particle """

        units = bond_length.unit
        dist = unit.Quantity(0.0,units)

        """ 'move_direction_list' tracks the Cartesian
        directions in which we have attempted a particle placement,
        where: x=0,y=1,z=2. """

        move_direction_list = []

        # Assign the parent coordinates as the initial coordinates for a trial particle

        trial_coordinates = parent_coordinates.__deepcopy__(memo={})
        ref = parent_coordinates.__deepcopy__(memo={})

        for direction in range(3):
            move_direction = random.randint(0,2)
            while move_direction in move_direction_list:
              move_direction = random.randint(0,2)

            if float(round(bond_length._value**2.0,4)-round(dist._value**2.0,4)) < 0.0:

              print("The bond length is: "+str(round(bond_length._value**2.0,4)))
              print("The distance is: "+str(round(dist._value**2.0,4)))
              print("The parent coordinates are: "+str(ref))
              print("The trial coordinates are: "+str(trial_coordinates))
              print("Error: new particles are not being assigned correctly.")
              exit()

            if direction == 2:
              trial_coordinates = get_move(trial_coordinates,move_direction,dist,bond_length,finish_bond=True)

            else:
              trial_coordinates = get_move(trial_coordinates,move_direction,dist,bond_length)

            move_direction_list.append(move_direction)

            print("Parent coordinates are: "+str(ref))
            print("Trial coordinates are: "+str(trial_coordinates))
            dist = distance(ref,trial_coordinates)
            print(direction)
            print(dist)

        if round(dist._value,4) < round(bond_length._value,4):

           print("Error: particles are being placed at a distance different from the bond length")
           print("Bond length is: "+str(bond_length))
           print("The particle distance is: "+str(dist))
           print(ref)
           print(trial_coordinates)
           exit()

        return(trial_coordinates)

def distances(interaction_list,positions):
        """
        Calculate the distances between a trial particle ('new_coordinates')
        and all existing particles ('existing_coordinates').

        Parameters
        ----------

        new_coordinates: Positions for a single trial particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        existing_coordinates: Positions for a single trial particle
        ( np.array( float * unit.angstrom ( shape = num_particles x 3 ) ) )

        Returns
        -------

        distances: List of the distances between all nonbonded particles.
        ( list ( float * simtk.unit.distance ( length = # nonbonded_interactions ) ) )

        """

        distance_list = []

#        print("In 'distances()' the positions are:"+str(positions))

        for interaction in interaction_list:
            if interaction[0] < len(positions) and interaction[1] < len(positions):
#             print("The distance between particles: "+str(interaction[0])+" and "+str(interaction[1])+" is "+str(distance(positions[interaction[0]],positions[interaction[1]])))
             distance_list.append(distance(positions[interaction[0]],positions[interaction[1]]))

        return(distance_list)

def collisions(distance_list,distance_cutoff):
        """
        Determine whether there are any collisions between non-bonded
        particles, where a "collision" is defined as a distance shorter
        than the user-provided 'bond_length'.

        Parameters
        ----------

        distances: List of the distances between all nonbonded particles.
        ( list ( float * simtk.unit.distance ( length = # nonbonded_interactions ) ) )

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        collision: Logical variable stating whether or not the model has
        bead collisions.
        default = False

        """

        collision = False

#        print("The nonbonded cutoff distance is: "+str(sigma._value))
        if len(distance_list) > 0:

          for distance in distance_list:

            if round(distance._value,4) < round(distance_cutoff._value,4):

              collision = True

        return(collision)

def assign_position_lattice_style(cgmodel,positions,distance_cutoff,bead_index,parent_index):
        """
        Assign random position for a bead

        Parameters
        ----------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """
        if bead_index == 0:
           print("This is the first bead.  Using the origin for coordinates.")
           success = True
           return(positions,success)

        units = cgmodel.bond_lengths['bb_bb_bond_length'].unit       
        if parent_index == -1:
               parent_index = bead_index - 1

        positions_before = positions.__deepcopy__(memo={})
#        print("Assigning coordinates for bead #"+str(bead_index))
#        print("(Indexing starts at: '1')")
#        print("Parent bead index is: "+str(parent_index-1))
#        print("Parent coordinates are: "+str(positions[parent_index-1]))
#        print("Child bead index is: "+str(bead_index-1))
#        print("Child coordinates are: "+str(positions[bead_index-1]))

#        print("Before assigning new coordinates the positions are: "+str(positions))
        bond_list = [[bond[0]-1,bond[1]-1] for bond in cgmodel.get_bond_list()]
       
#        print("The bond list is: "+str(bond_list))
        exclusion_list = []
        for i in range(cgmodel.num_beads):
          for j in range(i+1,cgmodel.num_beads):
            for bond_1 in range(len(bond_list)):
              if i in bond_list[bond_1]:
                if bond_list[bond_1][0] == i:
                  k = bond_list[bond_1][1]
                else:
                  k = bond_list[bond_1][0]
                for bond_2 in range(bond_1+1,len(bond_list)):
                  if j in bond_list[bond_2] and k in bond_list[bond_2]:
                    
                    if [i,j] not in exclusion_list and [j,i] not in exclusion_list:
                      if [i,j] not in bond_list and [j,i] not in bond_list:
                        exclusion_list.append([i,j])
#        print("The exclusion list is:"+str(exclusion_list))

        full_nonbonded_list = []
        for i in range(cgmodel.num_beads):
          for j in range(i+1,cgmodel.num_beads):
            if [i,j] not in bond_list and [j,i] not in bond_list:
              if [i,j] not in exclusion_list and [j,i] not in exclusion_list:
                full_nonbonded_list.append([i,j])
        
#        print("The list of non-bonded interactions is:"+str(full_nonbonded_list))
#        exit()

        parent_coordinates = positions[parent_index].__deepcopy__(memo={})
#        exit()
        new_coordinates = parent_coordinates.__deepcopy__(memo={})
        success = False
#        best_attempt = parent_coordinates.__deepcopy__(memo={})
#        best_distances = distances(full_nonbonded_list,positions[0:bead_index-1])
#        print(best_distances)
        move_direction_list = []
        while len(move_direction_list) < 6 and not success:
#           print(move_direction_list)
           bond_length = cgmodel.bond_lengths['bb_bb_bond_length']
#           bond_length = get_bond_length(cgmodel,parent_index,bead_index)
#           print("Before calling 'attempt_lattice_move' the value for new coordinates is:")
#           print(new_coordinates)
           new_coordinates,move_direction_list = attempt_lattice_move(parent_coordinates,bond_length,move_direction_list)
#           exit()
#           print("After calling 'attempt_lattice_move' the value for new coordinates is:")
#           print(new_coordinates)

           test_positions = positions[0:bead_index+1].__deepcopy__(memo={})
           test_positions[bead_index] = new_coordinates
#           print("Getting the non-bonded distances for these trial positions: ")
#           print(new_coordinates)
#           print(test_positions)
#           exit()
#           print(test_positions)
#           print("and this interaction list:")
#           print(full_nonbonded_list)
           distance_list = distances(full_nonbonded_list,test_positions)
#           print(full_nonbonded_list)
#           print("The distances are:"+str([distance._value for distance in distance_list]))
#           exit()
           if len(distance_list) == 0:
             success = True
             break
           if not collisions(distance_list,distance_cutoff) and len(distance_list) > 0: 
             success = True
#             print("The current list of nonbonded particle distances is: "+str(distance_list))
             break 

#           if collisions(distance_list,distance_cutoff) and len(distance
#             if min(distance_list) > min(best_distances):
#               best_attempt = test_positions.__deepcopy__(memo={})
#               best_distances = distance_list.__deepcopy__(memo={})
               
        if success:
          positions[bead_index] = new_coordinates
          if all(positions[bead_index] == positions_before[bead_index]):
            print("ERROR: the input coordinates are the same as the output coordinates in 'assign_position_lattice_style()'.")
            exit()
        return(positions,success)

def assign_position(positions,bond_length,sigma,bead_index,parent_index):
        """
        Assign random position for a bead

        Parameters
        ----------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """
        if bead_index == 1:
           success = True
           return(positions,success)

        units = bond_length.unit       
        if parent_index == -1:
               parent_index = bead_index - 1

        parent_coordinates = positions[parent_index-1]

        new_coordinates = unit.Quantity(np.zeros(3), units)
        success = False
        attempts = 0

        while not success:

           new_coordinates = attempt_move(parent_coordinates,bond_length)

           distance_list = distances(new_coordinates,positions)

           if not collisions(distance_list,sigma): success = True

           if attempts > 100:
            print([distance._value for distance in distance_list])
            return(positions,success)

           attempts = attempts + 1

        positions[bead_index-1] = new_coordinates

        return(positions,success)


def get_structure_from_library( cgmodel ):
        """
        Given a coarse grained model class object, this function retrieves
        a set of positions for the model from the ensemble library, in:
        '../foldamers/ensembles/${backbone_length}_${sidechain_length}_${sidechain_positions}'
        If this coarse grained model does not have an ensemble library, an 
        error message will be returned and we will attempt to assign 
        positions at random with 'random_positions()'.

        Parameters
        ----------

        cgmodel: CGModel() class object.

        Returns
        -------
        
        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """

        # In its current form this subroutine does not save the positions for random configurations we generate from heteropolymers.  It only saves the positions for homopolymers.

        if len(cgmodel.monomer_types) > 1:

          print("The foldamers ensemble library does not currently store conformations for polymers composed of more than one unique monomer.\n")
          print("The 'random_positions()' subroutine will be called instead, with 'use_library'=False.")
          positions = random_positions(cgmodel,use_library=False)
          return(positions)

        else:

          monomer_type = cgmodel.monomer_types[0]
          ensembles_directory = str(str(__file__.split('src/utilities/util.py')[0])+"ensembles")
          if not os.path.exists(ensembles_directory):
            os.mkdir(ensembles_directory)
          model_directory = str(str(ensembles_directory)+"/"+str(monomer_type['backbone_length'])+"_"+str(monomer_type['sidechain_length'])+"_"+str(monomer_type['sidechain_positions']))
          if not os.path.exists(model_directory):
            os.mkdir(model_directory)

          # We determine a suitable name for the ensemble directory by combining the 'bb_bb_bond_length', 'bb_sc_bond_length', and 'sc_sc_bond_length' into a single string:
          ens_str = [monomer_type['bond_lengths']['bb_bb_bond_length']._value,monomer_type['bond_lengths']['bb_sc_bond_length']._value,monomer_type['bond_lengths']['sc_sc_bond_length']._value]
          ensemble_directory = str(str(model_directory)+"/bonds_"+str(ens_str[0])+"_"+str(ens_str[1])+"_"+str(ens_str[2]))
          generate_ensemble = False
          if not os.path.exists(ensemble_directory):
            os.mkdir(ensemble_directory)
            generate_ensemble = True
#            positions = random_positions(cgmodel,use_library=False)
          pdb_list = []
          for file in os.listdir(ensemble_directory):
            if file.endswith('.pdb'):
              pdb_list.append(str(str(ensemble_directory)+"/"+str(file)))
          if len(pdb_list) < 100:
            generate_ensemble = True

          if generate_ensemble:
            print("The foldamers ensemble library only contains "+str(len(pdb_list))+" structures with these settings.\n")
            print("The 'random_positions()' subroutine will be called instead, with 'use_library'=False,")
            print("in order to generate a total of "+str(100)+" configurations for the database,")
            print("before a specific configuration is chosen to assign random positions for this model.")
            index = 1
            while index <= 100:
              file_name = str(ensemble_directory+"/cg"+str(index)+".pdb")
              if not os.path.exists(file_name):
                cgmodel.positions = random_positions(cgmodel,use_library=False)
               #print("Minimizing the structure.")
                positions_before = cgmodel.positions.__deepcopy__(memo={})
           
                positions_after,energy = minimize_structure(cgmodel.topology,cgmodel.system,cgmodel.positions)
                cgmodel.positions = positions_after.in_units_of(unit.angstrom)
                positions_after = cgmodel.positions.__deepcopy__(memo={})
                if all([all(positions_before[i] == positions_after[i]) for i in range(0,len(positions_before))]):
                  print("ERROR: these random positions were not suitable for an initial minimization attempt.")
                  print("NOTE: this routine will run continuously, unless the user interrupts with the keyboard.")
                  print("If we are attempting to build a file with the same name index, repeatedly, then there is probably")
                  print("something wrong with our model/parameter settings.")
                  continue

                else:

                  print("Adding a configuration to the ensemble library:")
                  print(file_name)
                  write_pdbfile_without_topology(cgmodel,file_name,energy=energy)
                  index = index + 1

              else:
                
                index = index + 1

        pdb_list = []
        for file in os.listdir(ensemble_directory):
          if file.endswith('.pdb'):
            pdb_list.append(str(str(ensemble_directory)+"/"+str(file)))
        random_file = pdb_list[random.randint(0,len(pdb_list)-1)]
        print("Using the positions found in:\n")
        print(str(random_file)+"\n")
        pdb_mm_obj = PDBFile(random_file)
        positions = pdb_mm_obj.getPositions()

        return(positions)

def random_positions( cgmodel,max_attempts=1000,use_library=True ):
        """
        Assign random positions for all beads in a coarse-grained polymer.

        Parameters
        ----------

        cgmodel: CGModel() class object.

        max_attempts: The maximum number of times that we will attempt to build
        a coarse grained model with the settings in 'cgmodel'.
        default = 1000

        use_library: A logical variable determining if we will generate a new
        random structure, or take a random structure from the library in the following path:
        '../foldamers/ensembles/${backbone_length}_${sidechain_length}_${sidechain_positions}'
        default = True
        ( NOTE: By default, if use_library = False, new structures will be added to the
          ensemble library for the relevant coarse grained model.  If that model does not
          have an ensemble library, one will be created. )

        Returns
        -------
        
        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """
        total_attempts = 0
        if use_library:
            print("Attempting to find a suitable random starting configuration in the foldamers structural database.\n")
            positions = get_structure_from_library(cgmodel)
#            print(positions)
#            print(cgmodel.positions)
#            if all([all(positions[i] == cgmodel.positions[i]) for i in range(0,len(positions))]): break
            return(positions)

        units = cgmodel.bond_lengths['bb_bb_bond_length'].unit
        positions = np.zeros([cgmodel.num_beads,3]) * units
        bond_list = [[bond[0]-1,bond[1]-1] for bond in cgmodel.get_bond_list()]
        total_attempts = 0
        distance_cutoff = 1.4 * cgmodel.bond_lengths['bb_bb_bond_length']

        lattice_style = True
        bond_index = 0
        stored_positions = positions[0].__deepcopy__(memo={})
        while total_attempts < max_attempts and len(stored_positions) != len(positions):
#         print("Beginning build attempt #"+str(total_attempts+1)+".")
         index = bond_index
         while index in range(bond_index,len(bond_list)):
          particle_index = bond_list[index][1]
#          print("using the bond between particles "+str(bond_list[index][0])+" and "+str(bond_list[index][1]))
#          print("to assign coordinates for particle "+str(particle_index))
          bond = bond_list[index]
          positions = np.zeros([cgmodel.num_beads,3]) * units
          shape = stored_positions.shape
          if shape == '(3,)':
            positions[0] = stored_positions.__deepcopy__(memo={})
          else:
            positions[0:len(stored_positions)] = stored_positions.__deepcopy__(memo={})
          if lattice_style:
#            print("BEFORE assigning new coordinates the positions are:"+str(positions))
#            print("Assigning positions using bond between particles:")
#            print(str(bond[0])+" and "+str(bond[1]))
            positions,placement = assign_position_lattice_style(cgmodel,positions,distance_cutoff,bond[1],bond[0])
#            print("AFTER assigning new coordinates the positions are:"+str(positions))
          else:
            positions,placement = assign_position(positions,cgmodel.bond_length,distance_cutoff,bond[1],bond[0])
          if not placement:
            if particle_index < 3:
              print("ERROR: having difficulty assigning random coordinates for this model,")
              print("even with a small number of particles ( <= 4 ).")
              print("The current positions are: "+str(positions))
              print("Please check the parameter values and other settings for your coarse grained model.")
              print("Then check the 'random_positions()' subroutine for errors.")
              exit()
            bond_index = round(random.triangular(2,index,index))
#            print("Restarting random configuration build from bond #"+str(bond_index)) 
            total_attempts = total_attempts + 1
            bond = bond_list[bond_index]
            particle_index = min(bond) + 1

            stored_positions = positions[0:particle_index].__deepcopy__(memo={})
#            print("The new positions are:"+str(stored_positions)+"\n")
            index = len(bond_list) + 1
          if placement:
            positions_before_appension = stored_positions.__deepcopy__(memo={})
            length_positions_before_appension = len(stored_positions)
            if particle_index > 0:
             stored_positions = positions[0:particle_index+1].__deepcopy__(memo={})
            else:
             stored_positions = np.zeros([1,3]) * units
             stored_positions[0] = positions[0].__deepcopy__(memo={})
            positions_after_appension = stored_positions
            length_positions_after_appension = len(stored_positions)
            shape = str(stored_positions.shape)
            if str(positions_after_appension.shape) not in ['(3,)','(1, 3)','(2, 3)']:
             if length_positions_before_appension != length_positions_after_appension - 1:
              print("ERROR: coordinates are not being added correctly in the 'random_positions()' subroutine.")
              print("The positions before appension are:"+str(positions_before_appension))
              print("The positions after appension are:"+str(positions_after_appension))
              exit()
            index = index + 1
#            bond_index = index + 1
          if len(stored_positions) == len(positions):
           
           break
#          print("The stored positions are:"+str(stored_positions))
        print("Successfully built a random structure with the following positions:\n")
        print(positions)
#        print("Double-checking that this structure has no particle collisions...\n")
        full_nonbonded_list = []
        for i in range(cgmodel.num_beads):
          for j in range(i+1,cgmodel.num_beads):
            if [i,j] not in bond_list:
              full_nonbonded_list.append([i,j])
        distance_list = distances(full_nonbonded_list,positions)
#           print(distance_list)
        if not collisions(distance_list,distance_cutoff):
#          print("No collisions found.")
          print("The shortest nonbonded particle distance is:\n")
          print(str(min(distance_list))+"\n")
          print("A nonbonded particle cutoff distance of "+str(distance_cutoff))
          print("was applied as criteria for successful structure generation.\n")
          return(positions)
        else:
          print("Error: A model was successfully built, however,")
          print("particle collisions were detected.\n")
          print("The shortest nonbonded particle distance is:")
          print(str(min(distance_list)))
          collision = full_nonbonded_list[distance_list.index(min(distance_list))]
          print("(between particles "+str(collision[0])+" and "+str(collision[1])+")")
          print("The nonbonded particle cutoff distance used for")
          print("random structure generation is set to:"+str(distance_cutoff))
#          exit()
          print("Going to attempt to generate another random structure.")
          print("This will continue until the user issues a disruption command with the keyboard. (Ctrl + c)")
          random_positions(cgmodel,max_attempts=1000,use_library=False)
        return(positions)

def distance(positions_1,positions_2):
        """
        Construct a matrix of the distances between all particles.

        Parameters
        ----------

        positions_1: Positions for a particle
        ( np.array( length = 3 ) )

        positions_2: Positions for a particle
        ( np.array( length = 3 ) )

        Returns
        -------

        distance
        ( float * unit )
        """

        direction_comp = np.zeros(3) * positions_1.unit

        for direction in range(len(direction_comp)):
          direction_comp[direction] = positions_1[direction].__sub__(positions_2[direction])

        direction_comb = np.zeros(3) * positions_1.unit.__pow__(2.0)
        for direction in range(3):
          direction_comb[direction] = direction_comp[direction].__pow__(2.0)

        sqrt_arg = direction_comb[0].__add__(direction_comb[1]).__add__(direction_comb[2])

        value = math.sqrt(sqrt_arg._value)
        units = sqrt_arg.unit.sqrt()
        distance = unit.Quantity(value=value,unit=units)

        return(distance)


def distance_matrix(positions):
        """
        Construct a matrix of the distances between all particles.

        Parameters
        ----------

        positions: Positions for an array of particles.
        ( np.array( num_particles x 3 ) )

        Returns
        -------

        distance_matrix: Matrix containing the distances between all beads.
        ( np.array( num_particles x 3 ) )
        """

        distance_matrix = np.array([[0.0 for index in range(0,len(positions))] for index in range(0,len(positions))])

        for index_1 in range(0,len(positions)):
          for index_2 in range(0,len(positions)):
            distance_matrix[index_1][index_2] = get_distance(positions[index_1],positions[index_2])

        return(distance_matrix)   
