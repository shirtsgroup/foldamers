
# =============================================================================================
# 1) PYTHON PACKAGE IMPORTS
# =============================================================================================

# System packages
import numpy as np
import math, random
from simtk import unit
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

def unit_sqrt(simtk_quantity):
        """
        Returns the square root of a simtk 'Quantity'.

        Parameters
        ----------

        simtk_quantity: A 'Quantity' object, as defined in simtk.
        ( float * unit )

        Returns
        -------

        answer: Square root of a simtk_quantity.

        """

        value = math.sqrt(simtk_quantity._value)
        answer = unit.Quantity(value,simtk_quantity.unit)

        return(answer)

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

        if positions != None:

           first_bead = False

        return(first_bead)

def single_bead(positions):
        """
        Determine if we have one particle in positions

        Parameters
        ----------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( float * unit ( shape = num_beads x 3 ) ) )

        Returns
        -------

        single_bead: Logical variable stating if this is the first particle.

        """
        single_bead = False
        if str(positions.shape) == '(3,)':
          single_bead = True

        return(single_bead)

def append_position(positions,new_coordinates):
        """
        Updates a set of input coordinates with 'new_coordinate' in the
        cartesian coordinate direction indexted by 'direction'.

        Parameters
        ----------

        new_coordinates: Cartesian coordinates for a particle
        ( np.array( float * unit ( length = 3 ) ) )

        direction: Cartesian direction index for particle placement, 
        where: x=0,y=1,z=2. 
        ( integer )

        trial_coordinates: Existing cartesian coordinates for the particle
        we are updating.
        ( np.array( float * unit ( length = 3 ) ) )
        Optional, default = None

        Returns
        -------

        trial_coordinates: Updated coordinates for the particle.

        """

        units = new_coordinates.unit
        if positions == None:

          new_positions = new_coordinates

        else:

          new_positions = unit.Quantity(np.zeros((2,3)),units)
          if single_bead(positions):
            new_positions = unit.Quantity(np.zeros((2,3)),units)
            new_positions[0] = positions
            new_positions[1] = new_coordinates

          if not single_bead(positions):
            new_positions = unit.Quantity(np.zeros((len(positions)+1,3)),units)
            for particle in range(len(positions)):
              new_positions[particle] = positions[particle]
            new_positions[len(new_positions)-1] = new_coordinates

        return(new_positions)

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
        max_step_size = unit_sqrt(bond_length.__pow__(2.0).__sub__(distance.__pow__(2.0)))

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
        trial_coordinates[move_direction] = step * trial_coordinates.unit
        return(trial_coordinates)

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

        trial_coordinates = np.array([parent_coordinates[i]._value for i in range(3)]) * parent_coordinates.unit
        ref = np.array([parent_coordinates[i]._value for i in range(3)]) * parent_coordinates.unit

        for direction in range(3):

            move_direction = random.randint(0,2)
            while move_direction in move_direction_list:
              move_direction = random.randint(0,2)

            if float(round(bond_length._value**2.0,4)-round(dist._value**2.0,4)) < 0.0:

              print("Error: new particles are not being assigned correctly.")
              exit()

            if direction == 2:
              trial_coordinates = get_move(trial_coordinates,move_direction,dist,bond_length,finish_bond=True)

            else:
              trial_coordinates = get_move(trial_coordinates,move_direction,dist,bond_length)

            move_direction_list.append(move_direction)

            dist = distance(ref,trial_coordinates)

        if round(dist._value,4) < round(bond_length._value,4):

           print("Error: particles are being placed at a distance different from the bond length")
           print("Bond length is: "+str(bond_length))
           print("The particle distance is: "+str(dist))
           print(ref)
           print(trial_coordinates)
           exit()

        return(trial_coordinates)

def non_bonded_distances(new_coordinates,existing_coordinates):
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

        distances = []

        if first_bead(existing_coordinates):

          return(distances)

        else:

          if single_bead(existing_coordinates):
              distances.append(distance(new_coordinates,existing_coordinates))

          else:

            for particle in range(0,len(existing_coordinates)):
              distances.append(distance(new_coordinates,existing_coordinates[particle]))

        return(distances)

def collisions(distances,bond_length):
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

        if len(distances) > 0:

          for distance in distances:

            if round(distance._value,4) < round(bond_length._value/2.0,4):

              collision = True

        return(collision)

def assign_position(positions,bond_length,bead_index,parent_index=-1):
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

        units = bond_length.unit       
        if parent_index == -1:
               parent_index = len(positions) - 1

        parent_coordinates = positions[parent_index]

        new_coordinates = unit.Quantity(np.zeros(3), units)
        success = False
        attempts = 0

        while not success:

           new_coordinates = attempt_move(parent_coordinates,bond_length)

           distances = non_bonded_distances(new_coordinates,positions)

           if not collisions(distances,bond_length): success = True

           if not success and attempts > 50:

            print("Error: maximum number of bead placement attempts exceeded")
            print("Double-check the bond lengths (and units) in your coarse-grained model.")
            return(positions,success)

           attempts = attempts + 1

        positions[bead_index] = new_coordinates
        success = True

        return(positions,success)

def random_positions( cgmodel ):
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

        positions = np.zeros([cgmodel.num_beads,3]) * cgmodel.bond_length.unit
        bond_list = cgmodel.get_bond_list()
        attempts = 0
        finished = False
        for bond in bond_list:
          while attempts < 100 and not finished:
           print("Assigning positions for particle # "+str(bond[0]))
           positions,placement = assign_position(positions,cgmodel.bond_length,bond[0],bond[1])
           attempts = attempts + 1
           if placement: 
          if bond == bond_list[len(bond_list)-1] and placement == True: finished = True

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
