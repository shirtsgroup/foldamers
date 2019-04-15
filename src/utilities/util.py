## No default python environment

# This script generates random coordinates for a CG polymer

# =============================================================================================
# 1) PYTHON PACKAGE IMPORTS
# =============================================================================================

# System packages
import numpy as np
import math, random
from simtk import unit
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

        sqrt: Square root of a simtk_quantity.

        """

        value = math.sqrt(simtk_quantity._value)

        return(value)

def append_position(positions,new_coordinate):
        """
        Updates a set of input coordinates with 'new_coordinate' in the
        cartesian coordinate direction indexted by 'direction'.

        Parameters
        ----------

        new_coordinate: Cartesian coordinates for a particle
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

        positions_values = positions._value
        units = new_coordinate.unit

        new_coordinate = np.array([float(coord) for coord in new_coordinate._value])
        positions_values = np.vstack([positions_values,new_coordinate])

        new_positions = unit.Quantity(positions_values,units)

        return(new_positions)

def update_trial_coordinates(move,trial_coordinates=None):
        """
        Updates 'trial_coordinates by adding the coordinates in 'move'.

        Parameters
        ----------

        move: Cartesian coordinates for a particle
        ( np.array( float * unit ( length = 3 ) ) )

        trial_coordinates: Existing cartesian coordinates for the particle
        we are updating.
        ( np.array( float * unit ( length = 3 ) ) )
        Optional, default = None

        Returns
        -------

        new_coordinates: Updated coordinates for the particle.

        """
        units = move.unit

        if trial_coordinates == None:

          trial_coordinates = unit.Quantity(np.zeros([3]), units)

        else:

          new_coordinates = unit.Quantity(np.zeros([3]), units)
          for direction in range(0,3):
            value = trial_coordinates[direction]._value
            new_coordinates[direction]._value = value

        return(new_coordinates)

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

        for value in positions._value:

          if type(value) == np.ndarray:

            first_bead = False

        return(first_bead)

def get_move(direction,step):
        """
        Given the cartesian coordinates for a particle ('move'),
        a 'step' (distance), and a 'direction' ( Index denoting
        x,y,z Cartesian direction), update the coordinates for
        the particle.

        Parameters
        ----------

        direction: Cartesian directions in which we have attempted 
        a particle placement, where: x=0,y=1,z=2. 
        ( integer )

        step: Number to add/subtract to the cartesian coordinates
        for direction 'direction' in 'move'
        ( float * simtk.unit.distance )

        Returns
        -------

        move: Updated positions for the particle
        ( np.array( float * unit.angstrom ( length = 3 ) ) )

        """

        move = unit.Quantity(np.zeros([3]),step.unit)

        if direction != 2:
          new_jump = random_sign(random.uniform(0.0,step._value))

        if direction == 2:
          new_jump = random_sign(step._value)

        values = [value for value in move._value]
        values[direction] = new_jump
        move = unit.Quantity(values,move.unit)
        return(move)

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

        units = parent_coordinates.unit
        dist = unit.Quantity(0.0,units)

        """ 'move_direction_list' tracks the Cartesian
        directions in which we have attempted a particle placement,
        where: x=0,y=1,z=2. """

        move_direction_list = []

        if first_bead(parent_coordinates):

          trial_coordinates = parent_coordinates

        if not first_bead(parent_coordinates):
          # Assign the parent coordinates as the initial coordinates for a trial particle
          trial_coordinates = unit.Quantity(parent_coordinates._value[len(parent_coordinates._value)-1],parent_coordinates.unit)

        ref = trial_coordinates

        for direction in range(0,3):

         move_direction = random.randint(0,2)
         while move_direction in move_direction_list:
          move_direction = random.randint(0,2)

         if float(round(bond_length._value**2.0,4)-round(dist._value**2.0,4)) < 0.0:

           print("Error: new particles are not being assigned correctly.")
           exit()

         step_arg = unit.Quantity(dist._value,units)
         step = unit.Quantity(unit_sqrt(step_arg),dist.unit)

         move = get_move(move_direction,step)
         move_direction_list.append(move_direction)

         trial_coordinates = update_trial_coordinates(move,trial_coordinates)

         dist = distance(ref,trial_coordinates)

         if round(dist._value,4) < round(bond_length._value,4):

           print("Error: particles are being placed at a distance different from the bond length")
           print("Bond length is: "+str(bond_length))
           print("The particle distance is: "+str(distance(trial_coordinates,ref)))
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

          for particle in range(0,len(existing_coordinates)):

            test_position = unit.Quantity(existing_coordinates._value[particle],existing_coordinates.unit)

            distances.append(distance(new_coordinates,test_position))

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

            if round(distance._value,4) < round(bond_length._value,4):

              collision = True

        return(collision)

def assign_position(positions,bond_length,parent_index=-1):
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

        units = positions.unit
        if len(positions) == 0:

          new_coordinates = unit.Quantity(np.zeros([3]), units)
          return
    
        else:

          if parent_index == -1:
            parent_coordinates = positions

          if parent_index == 0:
            parent_coordinates = positions

          if parent_index > 0:
            parent_coordinates = positions[parent_index-1]

        new_coordinates = unit.Quantity(np.zeros([3]), units)
        success = False
        attempts = 0

        while not success:

          new_coordinates = attempt_move(parent_coordinates,bond_length)

          distances = non_bonded_distances(new_coordinates,positions)

          if not collisions(distances,bond_length): success = True

          if not success and attempts > 1000000:

            print("Error: maximum number of bead placement attempts exceeded")
            print("Double-check the bond lengths (and units) in your coarse-grained model.")
            exit()

          attempts = attempts + 1

        positions = append_position(positions,new_coordinates)

        return(positions)

def assign_sidechain_beads( positions, sidechain_length, bond_length ):
        """
        Assign random position for all sidechain beads

        Parameters
        ----------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        sidechain_length: Number of beads in the sidechain
        portion of each (individual) monomer (integer), default = 1

        bond_length: Bond length for all beads that are bonded,
        ( float * simtk.unit.distance )
        default = 1.0 * unit.angstrom

        Returns
        -------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """

        for sidechain in range(0,sidechain_length):

          positions = assign_position(positions,bond_length)

        return(positions)

def assign_backbone_beads( positions, monomer_start, backbone_length, sidechain_length, sidechain_positions, bond_length ):
        """
        Assign random position for a backbone bead

        Parameters
        ----------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        monomer_start: Index of the bead to which we will bond this
        new backbone bead.
        ( integer )

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

        Returns
        -------

        positions: Positions for all beads in the coarse-grained model.
        ( np.array( num_beads x 3 ) )

        """

        for backbone_bead_index in range(0,backbone_length):

          if backbone_bead_index == 0:

            if not first_bead(positions):

              positions = assign_position(positions,bond_length,parent_index=monomer_start)

            else:

              positions = assign_position(positions,bond_length)

          # Assign side-chain beads

          if backbone_bead_index in sidechain_positions:

            positions = assign_sidechain_beads(positions,sidechain_length,bond_length)

        return(positions)


def random_positions( polymer_length, backbone_length, sidechain_length, sidechain_positions, bond_length, sigma ):
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

        positions = unit.Quantity(np.zeros([3]), unit.angstrom)

        for monomer in range(0,polymer_length):

          if monomer == 0:
            monomer_start = 0

          if monomer != 0:
            monomer_start = len(positions) - sidechain_length - 1

        # Assign backbone bead positions
        positions = assign_backbone_beads(positions,monomer_start,backbone_length,sidechain_length,sidechain_positions,bond_length)

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

        direction_comp = [0.0 for i in range(0,3)]

        for direction in range(len(direction_comp)):
          direction_comp[direction] = positions_1[direction].__sub__(positions_2[direction])

        direction_comb = [0.0 for i in range(0,3)]
        for direction in range(len(direction_comb)):
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
