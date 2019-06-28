from simtk import unit
import sys, os
from foldamers.src.utilities import util
from simtk import openmm as mm
from simtk.openmm.app.topology import Topology
from simtk.openmm.app.topology import Residue
import simtk.openmm.app.element as elem
from cg_openmm.src.build.cg_build import add_new_elements, build_topology, build_system
from itertools import chain, combinations, product

def get_parent_bead(cgmodel,monomer_index,bead_index,backbone_bead_index=None,sidechain_bead=False):
        """
        Determines the particle to which a given particle is bonded.  (Used for coarse grained model construction.)

        Parameters
        ----------

        cgmodel: CGModel() class object

        monomer_index: Index of the monomer the child particle belongs to.
        ( integer )
        Default = None

        bead_index: Index of the particle for which we would like to determine the parent particle it is bonded to.
        ( integer )
        Default = None

        backbone_bead_index: If this bead is a backbone bead, this index tells us its index (within a monomer) along the backbone
        ( integer )
        Default = None

        sidechain_bead: Logical variable stating whether or not this bead is in the sidechain.
        ( Logical )
        Default = False

        Returns
        -------

        parent_bead: Index for the particle that 'bead_index' is bonded to.
        ( Integer )

        """

        parent_bead = -1
        if bead_index != 1:
               if sidechain_bead == True:
#                print("Detecting type as sidechain when identifying parent.")
#                print(bead_index)
                parent_bead = bead_index - 1
#                print(parent_bead)
                return(parent_bead)
               else:
                 monomer_type = cgmodel.sequence[monomer_index]
                 if int(backbone_bead_index) in [monomer_type['sidechain_positions']]:
                  parent_bead = bead_index - monomer_type['sidechain_length'] - 1
                 else:
                  parent_bead = bead_index - 1
        if parent_bead == -1:
         print("ERROR: Not assigning particle indices correctly in get_parent_bead()")
         print("The bead index is: "+str(bead_index))
         exit()
        return(parent_bead)

def basic_cgmodel(polymer_length=8,backbone_length=1,sidechain_length=1,sidechain_positions=[0],mass=12.0 * unit.amu,bond_length=1.0 * unit.angstrom,sigma=2.5*unit.angstrom,epsilon=0.5 * unit.kilocalorie_per_mole,positions=None):

        """
        
        Parameters
        ----------

        :param polymer_length: Number of monomer units, default = 8
        :type polymer_length: integer

        :param backbone_length: Defines the number of beads in the backbone (assumes all monomers have the same backbone length), default = 1
        :type backbone_length: integer

        :param sidechain_length: Defines the number of beads in the sidechain ( assumes all monomers have the same sidechain length), default = 1
        :type sidechain_length: integer

        :param sidechain_positions: Defines the backbone bead indices upon which we will place the sidechains, default = [0]
        :type sidechain_positions: List( integer )

        :param mass: Mass for all coarse grained beads, default = 12.0 * unit.amu
        :type mass: float * simtk.unit

        :param bond_length: Defines the length for all bond types, default = 1.0 * unit.angstrom
        :type bond_length: float * simtk.unit
             
        :param sigma: Non-bonded bead Lennard-Jones equilibrium interaction distance, default = 2.5 * bond_length (for all particle interactions)
        :type sigma: float * simtk.unit

        :param epsilon: Non-bonded Lennard-Jones equilibrium interaction energy, default = 0.5 * unit.kilocalorie_per_mole
        :type espilon: float * simtk.unit

        :param positions: Positions for coarse grained particles in the model, default = None
        :type positions: np.array( float * simtk.unit ( shape = num_beads x 3 ) )

        Returns
        -------

        cgmodel: CGModel() class object

        """
        backbone_lengths = [backbone_length] # Number of backbone beads in unique monomer types
        sidechain_lengths = [sidechain_length] # Number of sidechain beads in unique monomer types
        masses = {'backbone_bead_masses': mass, 'sidechain_bead_masses': mass} # List of bead masses
        sigmas = {'bb_bb_sigma': sigma,'bb_sc_sigma': sigma,'sc_sc_sigma': sigma} # Lennard-Jones interaction distances.  List of unique interaction types
        bond_lengths = {'bb_bb_bond_length': bond_length,'bb_sc_bond_length': bond_length,'sc_sc_bond_length': bond_length} # bond length
        epsilons = {'bb_bb_eps': epsilon,'bb_sc_eps': epsilon,'sc_sc_eps': epsilon} # Lennard-Jones interaction strength.  List of unique interaction types
        cgmodel = CGModel(positions=positions,polymer_length=polymer_length,backbone_lengths=backbone_lengths, sidechain_lengths=sidechain_lengths, sidechain_positions = sidechain_positions, masses = masses, sigmas = sigmas, epsilons = epsilons, bond_lengths = bond_lengths)
        return(cgmodel)


class CGModel(object):
        """
        Returns
        -------

        cgmodel class object, detailing all of the properties for the coarse grained model.

        Parameters
        ----------

        :param positions: Positions for all of the particles, default = None
        :type positions: np.array( float * simtk.unit ( shape = num_beads x 3 ) )

        :param polymer_length: Length of the polymer, default = 8
        :type polymer_length: integer
      
        :param backbone_lengths: Defines the number of beads in the backbone for each monomer type, default = [1]
        :type backbone_lengths: List( integer )

        :param sidechain_lengths: Defines the number of beads in the sidechain for each monomer type, default = [1]
        :type sidechain_lengths: List( integer )

        :param sidechain_positions: Defines the backbone bead indices where sidechains are positioned, default = [0] (Place a sidechain on the first backbone bead in each monomer.)
        :type sidechain_positions: List( integer )

        :param masses: Masses of all particle types, default = 10.0 * unit.amu (for all particles)
        :type masses: dict( 'backbone_bead_masses': float * simtk.unit, 'sidechain_bead_masses': float * simtk.unit )

        :param sigmas: Non-bonded bead Lennard-Jones equilibrium interaction distances, default = 2.5 unit.angstrom (for all particles)
        :type sigmas: dict( 'bb_bb_sigma': float * simtk.unit,'bb_sc_sigma': float * simtk.unit,'sc_sc_sigma': float * simtk.unit} 

        :param epsilons: Non-bonded Lennard-Jones equilibrium interaction strengths, default = 0.5 * unit.kilocalorie_per_mole (for all particle interactions types)
        :type epsilons: dict( 'bb_bb_eps': float * simtk.unit,'bb_sc_eps': float * simtk.unit,'sc_sc_eps': float * simtk.unit )

        :param bond_lengths: Bond lengths for all bonds, default = 1.0 unit.angstrom
        :type bond_lengths: dict( 'bb_bb_bond_length': float * simtk.unit,'bb_sc_bond_length': float * simtk.unit,'sc_sc_bond_length': float * simtk.unit )

        :param bond_force_constants: Bond force constants for all bond types, default = 9.9e9 ( implied units are: kJ/mol/A^2 )
        :type bond_force_constants: dict( 'bb_bb_bond_k': float,'bb_sc_bond_k': float, 'sc_sc_bond_k': float )

        :param charges: Charges for all particles, default = 0.0 (for all particles)
        :type charges: dict( 'backbone_bead_charges': float * simtk.unit,'sidechain_bead_charges': float * simtk.unit )

        :param num_beads: Total number of particles in the coarse grained model, default = 16 (The total number of particles in a length=8 1-1 coarse-grained model)
        :type num_beads: integer

        :param system: OpenMM System() object, which stores the forces for the coarse grained model, default = None
        :type system: OpenMM System() class object

        :param topology: OpenMM Topology() object, which stores bonds, angles, and other structural attributes of the coarse grained model, default = None
        :type topology: OpenMM Topology() class object

        :param constrain_bonds: Logical variable determining whether bond constraints are applied during a simulation of the energy for the system, default = True
        :type constrain_bonds: Logical

        :param include_bond_forces: Include contributions from bond potentials when calculating the potential energy, default = True
        :type include_bond_forces: Logical

        :param include_nonbonded_forces: Include contributions from nonbonded interactions when calculating the potential energy, default = True
        :type include_nonbonded_forces: Logical

        :param include_bond_angle_forces: Include contributions from bond angle forces when calculating the potential energy, default = True
        :type include_bond_angle_forces: Logical

        :param include_torsion_forces: Include contributions from torsions when calculating the potential energy, default = True
        :type include_torsion_forces: Logical


        Attributes
        ----------

        **Attributes:**

        polymer_length : integer
                         Returns the number of monomers in the polymer/oligomer

        backbone_lengths : List( integers )
                           Returns a list of all unique backbone legnths (for individual monomers) in this model

        sidechain_lengths : List( integers )
                            Returns a list of all unique sidechain lengths (for individual monomers) in this model

        sidechain_positions : List( integers )
                              Returns a list of integers for all uniqye sidechain positions (along the backbone, for individual monomers) in this model

        masses : dict( float * simtk.unit )
                 Returns a list of the particle masses for all unique particle definitions in this model

        sigmas : dict( float * simtk.unit )
                 Returns a list of the Lennard-Jones nonbonded interaction distances for all unique particle interaction types

        epsilons : dict ( float * simtk.unit )
                   Returns a list of the Lennard-Jones nonbonded interaction strengths (well-depth) for all unique particle interaction types

        bond_lengths : dict ( float * simtk.unit )
                       Returns a list of the bond lengths for all unique bond definitions in the model

        nonbonded_interaction_list : List( List( integer, integer ) )
                                     Returns a list of the indices for particles that exhibit nonbonded interactions in this model

        bond_list : List( List( integer, integer ) )
                   Returns a list of paired indices for particles that are bonded in this model

        bond_angle_list : List( List( integer, integer, integer ) )
                          Returns a unique list of indices for all combinations of three particles that share a set of two bonds

        torsion_list: List( List( integer, integer, integer, integer ) )
                      Returns a unique list of indices for all ocmbinations of four particles that define a torsion (minimum requirement is that they share a set of three bonds)

        bond_force_constants : Dict( float )
                               Returns a dictionary with definitions for the bond force constants for all unique bond definitions

        bond_angle_force_constants: Dict( float )
                                    Returns a dictionary with definitions for the bond angle force constants for all unique bond angle definitions

        torsion_force_constants: Dict( float )
                                 Returns a dictionary with definitions for the torsion force constants for all unique torsion definitions

        equil_dihedral_angle : Dict( float )
                               Returns the equilibrium dihedral angle for all unique torsion definitions

        charges : Dict( float * simtk.unit )
                  Returns the charges for all unique particle definitions in this model

        num_beads : integer
                    Returns the number of particles in this model

        positions : np.array( float * simtk.unit ( shape = num_beads x 3 ) )
                    Returns the currently-stored positions for this model (if any)

        system : System() class object
                 Returns the currently-stored OpenMM System() object for this model (if any)

        topology : Topology() class object
                   Returns the currently-stored OpenMM Topology() object for this model (if any)

        constrain_bonds : Logical
                          Returns the current setting for bond constraints in the model

        include_bond_forces : Logical
                              Indicates if bond forces are currently included when calculating the energy

        include_nonbonded_forces : Logical
                                   Indicates if nonbonded interactions are currently included when calculating the energy

        include_bond_angle_forces : Logical
                                    Indicates if bond angle forces are currently included when calculating the energy

        include_torsion_forces : Logical
                                 Indicates if torsion potentials are currently included when calculating the energy
        
        """

        # Built in class attributes
        _BUILT_IN_REGIONS = ('polymer_length','backbone_lengths','sidechain_lengths','sidechain_positions','masses','sigmas','epsilons','bond_lengths','bond_force_constants','bond_angle_force_constants','torsion_force_constants','equil_dihedral_angle','charges','num_beads','positions','system','topology','constrain_bonds','bond_list','nonbonded_interaction_list','bond_angle_list','torsion_list','include_bond_forces','include_nonbonded_forces','include_bond_angle_forces','include_torsion_forces')

        def __init__(self,
                     positions = None,
                     polymer_length = 8,
                     backbone_lengths = [1],
                     sidechain_lengths = [1],
                     sidechain_positions = [0],
                     masses = {'backbone_bead_masses': 10.0 * unit.amu, 'sidechain_bead_masses': 10.0 * unit.amu}, 
                     sigmas = {'bb_bb_sigma': 2.5 * unit.angstrom,'bb_sc_sigma': 2.5 * unit.angstrom,'sc_sc_sigma': 2.5 * unit.angstrom},
                     epsilons = {'bb_bb_eps': 0.5 * unit.kilocalorie_per_mole,'bb_sc_eps': 0.5 * unit.kilocalorie_per_mole,'sc_sc_eps': 0.5 * unit.kilocalorie_per_mole}, 
                     bond_lengths = {'bb_bb_bond_length': 1.0 * unit.angstrom,'bb_sc_bond_length': 1.0 * unit.angstrom,'sc_sc_bond_length': 1.0 * unit.angstrom}, 
                     bond_force_constants = None, 
                     bond_angle_force_constants=None, 
                     torsion_force_constants=None, 
                     equil_bond_angle = None,
                     equil_dihedral_angle = None, 
                     charges = None, 
                     constrain_bonds = True,
                     include_bond_forces=True,
                     include_nonbonded_forces=True,
                     include_bond_angle_forces=True,
                     include_torsion_forces=True,
                     check_energy_conservation=True,
                     homopolymer=True):

          """
          Initialize variables that weren't provided
          """
          if bond_force_constants == None:
            bond_force_constants = {'bb_bb_bond_k': 9.9e9,'bb_sc_bond_k': 9.9e9, 'sc_sc_bond_k': 9.9e9}
          if bond_angle_force_constants == None:
            bond_angle_force_constants={'bb_bb_bb_angle_k': 200,'bb_bb_sc_angle_k': 200,'bb_sc_sc_angle_k': 200,'sc_sc_sc_angle_k': 200}
          if torsion_force_constants == None:
            torsion_force_constants={'bb_bb_bb_bb_torsion_k': 200,'bb_bb_bb_sc_torsion_k': 200,'bb_bb_sc_sc_torsion_k': 200, 'bb_sc_sc_sc_torsion_k': 200, 'sc_bb_bb_sc_torsion_k': 200, 'bb_sc_sc_bb_torsion_k': 200, 'sc_sc_sc_sc_torsion_k': 200}
          if equil_bond_angle == None:
            equil_bond_angle = 120
          if equil_dihedral_angle == None:
            equil_dihedral_angle = 180
          if charges == None:
            charges = {'backbone_bead_charges': 0.0 * unit.elementary_charge,'sidechain_bead_charges': 0.0 * unit.elementary_charge}

          """
          Initialize variables that were passed as input
          """

          self.polymer_length = polymer_length
          self.backbone_lengths = backbone_lengths
          self.sidechain_lengths = sidechain_lengths
          self.sidechain_positions = sidechain_positions
          self.bond_lengths = bond_lengths
          self.monomer_types = self.get_monomer_types()

          sequence = []
          if homopolymer == True:
           monomer_type = self.monomer_types[0]
           for monomer in range(self.polymer_length):
            sequence.append(monomer_type)
          else:
           for monomer in range(self.polymer_length):
            monomer_type_index = random.randint(0,len(self.monomer_types))
            sequence.append(monomer_types[monomer_type_index])

          self.sequence = sequence
          self.num_beads = self.get_num_beads()
          self.particle_list = self.get_particle_list()
          self.masses = masses
          self.sigmas = sigmas
          self.epsilons = epsilons
          self.bond_force_constants = bond_force_constants
          self.bond_angle_force_constants = bond_angle_force_constants
          self.equil_bond_angle = equil_bond_angle
          self.torsion_force_constants = torsion_force_constants
          self.equil_dihedral_angle = equil_dihedral_angle
          self.charges = charges

          if len(sigmas) != 0: include_nonbonded_forces = True
          if len(bond_force_constants) != 0: include_bond_forces = True
          if len(bond_angle_force_constants) != 0 and include_bond_angle_forces != False: include_bond_angle_forces = True
          if len(torsion_force_constants) != 0: include_torsion_forces = True

          self.include_bond_forces = include_bond_forces
          self.include_bond_angle_forces = include_bond_angle_forces
          self.include_nonbonded_forces = include_nonbonded_forces
          self.include_torsion_forces = include_torsion_forces
          self.check_energy_conservation = check_energy_conservation

          """
          Get bond, angle, and torsion lists.
          """
          self.bond_list = self.get_bond_list()
          self.nonbonded_interaction_list = self.get_nonbonded_interaction_list()
          self.bond_angle_list = self.get_bond_angle_list()
          self.torsion_list = self.get_torsion_list()
          self.constrain_bonds = constrain_bonds

          """
          Make a list of coarse grained particle masses:
          """
          list_of_masses = self.get_all_particle_masses()

          """
          Initialize new (coarse grained) particle types:
          """
          self.particle_types = add_new_elements(self,list_of_masses)

          self.system = build_system(self)
          
          self.topology = build_topology(self)

          if positions == None: self.positions = util.random_positions(self,use_library=True) 
          else: self.positions = positions

        def get_monomer_types(self):
          """
          Get a list of monomer dictionary objects for each unique monomer type.
          """
          monomer_name_modifier = ['A','B','C','D','E','F','G','H']
          monomer_types = []
          monomer_type_index = 0
          for backbone_length,sidechain_length in zip(self.backbone_lengths,self.sidechain_lengths):
            num_beads = backbone_length
            for sidechain_position in self.sidechain_positions:
             num_beads = num_beads + sidechain_length
            monomer_name = str('CG'+str(backbone_length)+str(sidechain_length))
            if monomer_name in monomer_types:
             modifier_index = 0
             while monomer_name in monomer_types:
              monomer_name = str('CG'+str(backbone_length)+str(sidechain_length)+str(modifier_index))
              modifier_index = modifier_index + 1
            monomer_type = {'monomer_name': monomer_name, 'backbone_length': backbone_length, 'sidechain_length': sidechain_length, 'sidechain_positions': sidechain_position, 'num_beads': num_beads, 'bond_lengths': self.bond_lengths}
            monomer_types.append(monomer_type)
          return(monomer_types)

        def get_num_beads(self):
          """
          Calculate the number of beads in our coarse grained model(s)
          """
          num_beads = 0
          for monomer_type in self.sequence:
           num_beads = num_beads + monomer_type['num_beads']
          return(num_beads)

        def get_particle_list(self):
          """
          Get a list of particles, where the indices correspond to those used in our system/topology
          """
          particle_list = []
          for monomer_type in self.sequence:
           cg_particle_index = 1
           for backbone_bead in range(monomer_type['backbone_length']):
            particle_symbol = str("B"+str(cg_particle_index))
            particle_list.append(particle_symbol)
            cg_particle_index = cg_particle_index + 1
            if type(monomer_type['sidechain_positions']) == int:
             sidechain_positions = [monomer_type['sidechain_positions']]
            else:
             sidechain_positions = monomer_type['sidechain_positions']
            if backbone_bead in sidechain_positions:
             for sidechain in range(monomer_type['sidechain_length']):
              particle_symbol = str("S"+str(cg_particle_index))
              particle_list.append(particle_symbol)
              cg_particle_index = cg_particle_index + 1
          return(particle_list)


        def get_bond_list(self):
          """
          Construct a bond list for the coarse grained model
          """
          bond_list = []
          bead_index = 1
          for monomer in range(len(self.sequence)):
            monomer_type = self.sequence[monomer]
            for backbone_bead in range(1,monomer_type['backbone_length']+1):

             if bead_index != 1:
#                             get_parent_bead(cgmodel,monomer_index,bead_index,backbone_bead_index=None,sidechain_bead=False):
#              print(bead_index)
              parent_index = get_parent_bead(self,monomer,bead_index,backbone_bead-1,sidechain_bead=False)
              if parent_index < 0 or bead_index < 0:
               print("Error: identifying parent index incorrectly when assigning a bond between backbone particles.")
               print("The bead index is: "+str(bead_index))
               print("The parent index is: "+str(parent_index))
               print("The backbone index is: "+str(backbone_bead-1))
               exit() 
              if parent_index < bead_index:
               bond_list.append([parent_index,bead_index])
              else:
               bond_list.append([bead_index,parent_index])
             bead_index = bead_index + 1
             
             if int(backbone_bead-1) in [monomer_type['sidechain_positions']]:
                for sidechain_bead in range(monomer_type['sidechain_length']):
                  parent_index = get_parent_bead(self,monomer,bead_index,backbone_bead-1,sidechain_bead=True)
                  if parent_index < 0 or bead_index < 0:
                   print("Error: identifying parent index incorrectly when assigning a bond with sidechain particles.")
                   print("The bead index is: "+str(bead_index))
                   print("The parent index is: "+str(parent_index))
                   print("The backbone index is: "+str(backbone_bead-1))
                   exit()                  
                  if parent_index < bead_index:
                   bond_list.append([parent_index,bead_index])
                  else:
                   bond_list.append([bead_index,parent_index])
                  bead_index = bead_index + 1
#          print(bond_list)
#          exit()
          return(bond_list)

        def get_nonbonded_interaction_list(self):
          """
          Construct a nonbonded interaction list for our coarse grained model
          """

          interaction_list = []
          bond_list = [[bond[0]-1,bond[1]-1] for bond in self.get_bond_list()]
          for particle_1 in range(self.num_beads):
               for particle_2 in range(self.num_beads):
                 if particle_1 != particle_2 and abs(particle_1 - particle_2) >= 3:
                   if [particle_1,particle_2] not in bond_list and [particle_2,particle_1] not in bond_list:
                     if [particle_1,particle_2] not in interaction_list:
                       if [particle_2,particle_1] not in interaction_list:
                         interaction_list.append([particle_1,particle_2])
                     if [particle_2,particle_1] not in interaction_list:
                       if [particle_1,particle_2] not in interaction_list:
                         interaction_list.append([particle_2,particle_1])
          return(interaction_list)


        def get_bond_angle_list(self):
          """
          Construct a list of indices for particles that define bond angles in our coarse grained model
          """

          bond_list = self.bond_list
          bond_angles = []
          for bond_1 in bond_list:
            bond_angle = [bond_1[0],bond_1[1]]
            for bond_2 in bond_list:
             if bond_2 != bond_1 and [bond_2[1],bond_2[0]] != bond_1:
              if bond_1[0] in bond_2 or bond_1[1] in bond_2:
               if bond_2[0] not in bond_angle:
                bond_angle.append(bond_2[0])
               if bond_2[1] not in bond_angle:
                bond_angle.append(bond_2[1])
             if len(bond_angle) == 3:
                 unique = True
                 for existing_bond_angle in bond_angles:
                  if all(bond_angle) in existing_bond_angle:
                   unique = False
                 if unique:
                   bond_angles.append(bond_angle)
                 bond_angle = bond_1

          return(bond_angles)


        def get_torsion_list(self):
          """
          Construct a torsion list for our coarse grained model
          """

          bond_list = self.bond_list
          torsions = []
          for bond_1 in bond_list:
            torsion = [bond_1[0],bond_1[1]]
            for bond_2 in bond_list:
             if bond_2 != bond_1 and [bond_2[1],bond_2[0]] != bond_1:
              if bond_1[0] in bond_2 or bond_1[1] in bond_2:
               if bond_2[0] not in torsion:
                torsion.append(bond_2[0])
               if bond_2[1] not in torsion:
                torsion.append(bond_2[1])

             for bond_3 in bond_list:
              if bond_3 != bond_1 and [bond_3[1],bond_3[0]] != bond_1:
                if bond_3 != bond_2 and [bond_3[1],bond_3[0]] != bond_2:
                  if bond_3[0] in bond_2 or bond_1[1] in bond_2:
                    if bond_2[0] not in torsion:
                      torsion.append(bond_2[0])
                    if bond_2[1] not in torsion:
                      torsion.append(bond_2[1])

             if len(torsion) == 3:
                 unique = True
                 for existing_torsion in torsions:
                  if all(torsion) in existing_torsion:
                   unique = False
                 if unique:
                   torsions.append(torsion)

          return(torsions)

        def get_particle_type(self,particle_index,particle_name=None):
          """
          Returns the name of a particle, given its index within the model

          Parameters
          ----------

          self: CGModel() class object

          particle_index: Index of the particle for which we would like to determine the type
          Type: int()

          Returns
          -------

          particle_type: 'backbone' or 'sidechain'
          Type: str()

          """
          if particle_name == None: particle_name = self.particle_list[particle_index]
          if 'B' in particle_name: particle_type = 'backbone'
          if 'S' in particle_name: particle_type = 'sidechain'

          return(particle_type)

        def get_particle_mass(self,particle_index):
          """
          Returns the mass for a particle, given its index.

          Parameters
          ----------

          self: CGModel() class object

          Returns
          -------

          Mass

          """
          particle_type = self.get_particle_type(particle_index)
          if particle_type == 'backbone': particle_mass = self.masses['backbone_bead_masses']
          if particle_type == 'sidechain': particle_mass = self.masses['sidechain_bead_masses']
          return(particle_mass)

        def get_particle_charge(self,particle_index):
          """
          Returns the charge for a particle, given its index.

          Parameters
          ----------

          self: CGModel() class object

          Returns
          -------

          Charge

          """
          particle_type = self.get_particle_type(particle_index)
          if particle_type == 'backbone': particle_charge = self.charges['backbone_bead_charges']
          if particle_type == 'sidechain': particle_charge = self.charges['sidechain_bead_charges']
          return(particle_charge)

        def get_sigma(self,particle_index,particle_type=None):
          """
          Returns the sigma value for a particle, given its index within the coarse grained model.

          Parameters
          ----------

          self: CGModel() class object

          Returns
          -------

          Sigma

          """

          if particle_type == None: particle_type = self.get_particle_type(particle_index)
          if particle_type == 'backbone': sigma = self.sigmas['bb_bb_sigma']
          if particle_type == 'sidechain': sigma = self.sigmas['sc_sc_sigma']
          return(sigma)

        def get_epsilon(self,particle_index,particle_type=None):
          """
          Returns the epsilon value for a particle, given its index.

          Parameters
          ----------

          self: CGModel() class object

          Returns
          -------

          Epsilon

          """
          if particle_type == None: particle_type = self.get_particle_type(particle_index)
          if particle_type == 'backbone': epsilon = self.epsilons['bb_bb_eps']
          if particle_type == 'sidechain': epsilon = self.epsilons['sc_sc_eps']
          return(epsilon)

        def get_all_particle_masses(self):
          """
          Returns a list of unique particle masses

          Parameters
          ----------

          self: CGModel() class object

          Returns
          -------

          List( unique particle masses )

          """
          list_of_masses = []
          list_of_masses.append(self.masses['backbone_bead_masses'])
          list_of_masses.append(self.masses['sidechain_bead_masses'])
          return(list_of_masses)

        def get_bond_length_from_names(self,particle_1_name,particle_2_name):
          """
          Determines the correct bond length for two particles, given their symbols.

          Parameters
          ----------

          cgmodel: CGModel() class object

          particle_1_name: Symbol for the first particle in the bond
          ( string )
          Default = None

          particle_2_name: Symbol for the second particle in the bond
          ( string )
          Default = None

          Returns
          -------

          bond_length: Bond length for the bond defined by these two particles.
          ( simtk.unit.Quantity() )

          """
          if 'B' in particle_1_name: particle_1_type = 'backbone'
          else: particle_1_type = 'sidechain'
          if particle_1_type == 'backbone' and particle_2_type == 'backbone':
           bond_length = self.bond_lengths['bb_bb_bond_length']

          if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
           bond_length = self.bond_lengths['bb_sc_bond_length']

          if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
           bond_length = self.bond_lengths['bb_bb_bond_length']

          return(bond_length)

        def get_bond_length(self,particle_1_index,particle_2_index):
          """
          Determines the correct bond force constant for two particles

          Parameters
          ----------

          self: CGModel() class object

          particle_1_index: Index of the first particle in the bond
          ( integer )
          Default = None

          particle_2_index: Index of the second particle in the bond
          ( integer )
          Default = None

          Returns
          -------

          bond_length: Bond length for the bond defined by these two particles.
          ( simtk.unit.Quantity() )

          """
          if 'B' in self.particle_list[particle_1_index]: particle_1_type = 'backbone'
          else: particle_1_type = 'sidechain'

          if 'B' in self.particle_list[particle_2_index]: particle_2_type = 'backbone'
          else: particle_2_type = 'sidechain'

          if particle_1_type == 'backbone' and particle_2_type == 'backbone':
           bond_length = self.bond_lengths['bb_bb_bond_length']

          if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
           bond_length = self.bond_lengths['bb_sc_bond_length']

          if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
           bond_length = self.bond_lengths['bb_bb_bond_length']

          return(bond_length)

        def get_bond_force_constant(self,particle_1_index,particle_2_index):
          """
          Determines the correct bond force constant for two particles

          Parameters
          ----------

          cgmodel: CGModel() class object

          particle_1_index: Index of the first particle in the bond, default = None

          particle_2_index: Index of the second particle in the bond, default = None

          Returns
          -------

          bond_force_constant: Bond force constant for the bond defined by these two particles

          """
          if 'B' in self.particle_list[particle_1_index]: particle_1_type = 'backbone'
          else: particle_1_type = 'sidechain'

          if 'B' in self.particle_list[particle_2_index]: particle_2_type = 'backbone'
          else: particle_2_type = 'sidechain'

          if particle_1_type == 'backbone' and particle_2_type == 'backbone':
           bond_force_constant = self.bond_force_constants['bb_bb_bond_k']

          if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
           bond_force_constant = self.bond_force_constants['bb_sc_bond_k']

          if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
           bond_force_constant = self.bond_force_constants['bb_bb_bond_k']

          return(bond_force_constant)

        def get_bond_angle(self,particle_1_index,particle_2_index,particle_3_index):
          """
          Determines the correct equilibrium bond angle between three particles

          Parameters
          ----------

          self: CGModel() class object

          particle_1_index: Index of the first particle in the bond, default = None

          particle_2_index: Index of the second particle in the bond angle, default = None

          particle_3_index: Index of the third particle in the bond angle, default = None

          Returns
          -------

          bond_angle: Bond angle for the two bonds defined by these three particles.
          """

          bond_angle = cgmodel.equil_bond_angle

          return(bond_angle)

        def get_bond_angle_force_constant(self,particle_1_index,particle_2_index,particle_3_index):
          """
          Determines the correct equilibrium bond angle between three particles

          Parameters
          ----------

          self: CGModel() class object

          particle_1_index: Index of the first particle in the bond, default = None

          particle_2_index: Index of the second particle in the bond angle, default = None

          particle_3_index: Index of the third particle in the bond angle, default = None

          Returns
          -------

          bond_angle: Bond angle for the two bonds defined by these three particles.
          """
          particle_1_type = self.get_particle_type(particle_1_index)
          particle_2_type = self.get_particle_type(particle_2_index)
          particle_3_type = self.get_particle_type(particle_3_index)

          bond_angle_force_constant = 200

          return(bond_angle_force_constant)

        def get_torsion_force_constant(self,torsion):
          """         
          Determines the torsion force constant given a list of particle indices

          Parameters
          ----------

          cgmodel: CGModel() class object

          torsion: Indices of the particles in the torsion
          ( integer )
          Default = None

          Returns
          -------

          torsion_force_constant: Force constant for the torsion defined by the input particles.
          ( Integer )

          """
          particle_types = ['','','','']
          if 'B' in self.particle_list[torsion[0]]: particle_types[0] = 'backbone'
          else: particle_types[0] = 'sidechain'

          if 'B' in self.particle_list[torsion[1]]: particle_types[1] = 'backbone'
          else: particle_types[1] = 'sidechain'

          if 'B' in self.particle_list[torsion[2]]: particle_types[2] = 'backbone'
          else: particle_types[2] = 'sidechain'

          if 'B' in self.particle_list[torsion[3]]: particle_types[3] = 'backbone'
          else: particle_types[3] = 'sidechain'

          if particle_types[0] == 'sidechain':
           if particle_types[1] == 'backbone':
            if particle_types[2] == 'backbone':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['sc_bb_bb_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['sc_bb_bb_sc_torsion_k']
            if particle_types[2] == 'sidechain':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['bb_sc_sc_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['sc_bb_sc_sc_torsion_k']
           if particle_types[1] == 'sidechain':
            if particle_types[2] == 'backbone':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['sc_sc_bb_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['sc_sc_bb_sc_torsion_k']
            if particle_types[2] == 'sidechain':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['sc_sc_sc_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['sc_sc_sc_sc_torsion_k']
          if particle_types[0] == 'backbone':
           if particle_types[1] == 'backbone':
            if particle_types[2] == 'backbone':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['bb_bb_bb_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['bb_bb_bb_sc_torsion_k']
            if particle_types[2] == 'sidechain':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['bb_bb_sc_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['bb_bb_sc_sc_torsion_k']
           if particle_types[1] == 'sidechain':
            if particle_types[2] == 'backbone':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['bb_sc_bb_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['bb_sc_bb_sc_torsion_k']
            if particle_types[2] == 'sidechain':
             if particle_types[3] == 'backbone':
              torsion_force_constant = self.torsion_force_constants['bb_sc_sc_bb_torsion_k']
             if particle_types[3] == 'sidechain':
              torsion_force_constant = self.torsion_force_constants['bb_sc_sc_sc_torsion_k']
          return(torsion_force_constant)

