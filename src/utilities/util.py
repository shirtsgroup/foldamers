
# =============================================================================================
# 1) PYTHON PACKAGE IMPORTS
# =============================================================================================

# System packages
import numpy as np
import math, random
from simtk import unit
#from foldamers.src.cg_model.cgmodel import get_bond_length
from foldamers.src.cg_model.cgmodel import *

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

        for interaction in interaction_list:
            if interaction[0] < len(positions) and interaction[1] < len(positions):
             print("The distance between particles: "+str(interaction[0])+" and "+str(interaction[1])+" is "+str(distance(positions[interaction[0]],positions[interaction[1]])))
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
        if bead_index == 1:
           success = True
           return(positions,success)

        units = cgmodel.bond_lengths['bb_bb_bond_length'].unit       
        if parent_index == -1:
               parent_index = bead_index - 1

#        print("Parent bead index is: "+str(parent_index-1))
#        print("Parent coordinates are: "+str(positions[parent_index-1]))
#        print("Child bead index is: "+str(bead_index-1))
#        print("Child coordinates are: "+str(positions[bead_index-1]))
#        print("Positions are: "+str(positions))
        bond_list = [[bond[0]-1,bond[1]-1] for bond in cgmodel.get_bond_list()]
       
#        print("The bond list is: "+str(bond_list))
        full_nonbonded_list = []
        for i in range(cgmodel.num_beads):
         for j in range(i+1,cgmodel.num_beads):
           if [i,j] not in bond_list:
             full_nonbonded_list.append([i,j])
#        print(full_nonbonded_list)
        parent_coordinates = positions[parent_index-1].__deepcopy__(memo={})

        new_coordinates = unit.Quantity(np.zeros(3), units)
        success = False
        best_attempt = positions[parent_index-1].__deepcopy__(memo={})
        best_distances = distances(full_nonbonded_list,positions[0:bead_index-1])
#        print(best_distances)
        move_direction_list = []
        while len(move_direction_list) < 6 and not success:
           bond_length = cgmodel.bond_lengths['bb_bb_bond_length']
#           bond_length = get_bond_length(cgmodel,parent_index,bead_index)
           new_coordinates,move_direction_list = attempt_lattice_move(parent_coordinates,bond_length,move_direction_list)

           test_positions = positions.__deepcopy__(memo={})
           test_positions[bead_index-1] = new_coordinates
           print("Getting the non-bonded distances for these trial positions: ")
           print(test_positions[bead_index-1])
           distance_list = distances(full_nonbonded_list,test_positions[0:bead_index-1])
#           print(distance_list)
           if not collisions(distance_list,distance_cutoff): 
             success = True
             print("No collisions found")

           if collisions(distance_list,distance_cutoff) and len(distance_list) > 0:
             if min(distance_list) > min(best_distances):
               best_attempt = test_positions.__deepcopy__(memo={})
               best_distances = distance_list.__deepcopy__(memo={})
               
        if len(move_direction_list) == 6 and not success:
          print("Couldn't find a place to put a coarse grained bead.")
          exit()
        print("The new coordinates are: "+str(test_positions))
        positions = test_positions
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

def random_positions( cgmodel,max_attempts=1000 ):
        """
        Assign random positions for all beads in a coarse-grained polymer.

        Parameters
        ----------

        polymer_length: Number of monomer units (integer), default = 8
      
        backbone_length: Number of beads in the backbone 
        portion of each (individual) monomer (integer), default = 1

        sidechain_length: Number of beads in the sidechain
        portion of each (individual) monomer (integer), default = 1

        sidechain_positions: List of integers defining the backbone
        bead indices upon which we will place the sidechains,
        default = [0] (Place a sidechain on the backbone bead with
        index "0" (first backbone bead) in each (individual) monomer

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        sigma: Non-bonded bead Lennard-Jones interaction distances,
        ( float * simtk.unit.distance )
        default = 8.4 * unit.angstrom

        bb_bond_length: Bond length for all bonded backbone beads,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        bs_bond_length: Bond length for all backbone-sidechain bonds,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        ss_bond_length: Bond length for all beads within a sidechain,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------
        
        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """

        units = cgmodel.bond_lengths['bb_bb_bond_length'].unit
        positions = np.zeros([cgmodel.num_beads,3]) * units
        bond_list = cgmodel.get_bond_list()
        bond_index = 0
        total_attempts = 0
        distance_cutoff = 1.2 * cgmodel.bond_lengths['bb_bb_bond_length']

        lattice_style = True

        while bond_index in range(len(bond_list)) and total_attempts < max_attempts:
          bond = bond_list[bond_index]
          attempts = 0
          placement = False
          while attempts < max_attempts and not placement:
           if lattice_style:
            stored_positions = positions.__deepcopy__(memo={})
            positions,placement = assign_position_lattice_style(cgmodel,positions,distance_cutoff,bond[1],bond[0])
           else:
            positions,placement = assign_position(positions,cgmodel.bond_length,distance_cutoff,bond[1],bond[0])
           attempts = attempts + 1
          if not placement: 
           bond_index = bond_index - round(random.triangular(0,bond_index))
           total_attempts = total_attempts + 1
           positions = stored_positions
          if placement: bond_index = bond_index + 1
        if total_attempts >= max_attempts:
         print("Failed to build a CG model with "+str(cgmodel.backbone_length)+" backbone beads")
         print(" and "+str(cgmodel.sidechain_length)+" sidechain beads")

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
