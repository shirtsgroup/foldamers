from simtk import unit
import sys, os
from foldamers.src.utilities import util
from simtk import openmm as mm
import simtk.openmm.app.element as elem

def get_particle_masses(cgmodel):
        """
        Returns a list of unique particle masses

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        List( unique particle masses )

        """
        list_of_masses = []
        for backbone_bead in range(cgmodel.backbone_length):
            list_of_masses.append(cgmodel.mass)
            if backbone_bead in cgmodel.sidechain_positions:
              for sidechain in range(cgmodel.sidechain_length):
                 list_of_masses.append(cgmodel.mass)
        return(list_of_masses)

def add_new_elements(cgmodel,list_of_masses):
        """
        Adds new coarse grained particle types to OpenMM

        Parameters
        ----------

        cgmodel: CGModel() class object

        list_of_masses: List of masses for the particles we want to add to OpenMM

        """
        element_index = 117
        mass_index = 0
        cg_particle_index = 1
        for backbone_bead in range(cgmodel.backbone_length):

         particle_name = str("bb-"+str(cg_particle_index))
         particle_symbol = str("B"+str(cg_particle_index))
         if particle_symbol not in elem.Element._elements_by_symbol:
          elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
          element_index = element_index + 1
          cg_particle_index = cg_particle_index + 1
          mass_index = mass_index + 1
         if backbone_bead in cgmodel.sidechain_positions:
           for sidechain in range(cgmodel.sidechain_length):
            if particle_symbol not in elem.Element._elements_by_symbol:
             particle_name = str("sc-"+str(cg_particle_index))
             particle_symbol = str("S"+str(cg_particle_index))
             elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
             element_index = element_index + 1
             cg_particle_index = cg_particle_index + 1
             mass_index = mass_index + 1

        return

def get_parent_bead(cgmodel,bead_index,backbone_bead_index=None,sidechain_bead=False):
        """
        Determines the particle to which a given particle is bonded.  (Used for coarse grained model construction.)

        Parameters
        ----------

        cgmodel: CGModel() class object

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
                parent_bead = bead_index - 1
                return(parent_bead)
               else:
                 if backbone_bead_index - 1 in cgmodel.sidechain_positions:
                  parent_bead = bead_index - cgmodel.sidechain_length - 1
                 else:
                  parent_bead = bead_index - 1
        return(parent_bead)

def build_system(cgmodel):
        """
        Builds an OpenMM System() class object, given a CGModel() class object as input.

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        system: OpenMM System() class object

        """
        sigma = cgmodel.sigma.in_units_of(unit.nanometer)._value
        charge = cgmodel.charge._value
        epsilon = cgmodel.epsilon.in_units_of(unit.kilojoule_per_mole)._value
        bond_length = cgmodel.bond_length.in_units_of(unit.nanometer)._value

        # Create system
        system = mm.System()

        if cgmodel.include_bond_forces:
         # Create bond (harmonic) potentials
         bond_list = cgmodel.get_bond_list()
         new_bond_list = []
         bead_index = 1
         bond_force = mm.HarmonicBondForce()
         for bond in bond_list:
              new_bond = [bond[0]-1,bond[1]-1]
              new_bond_list.append(new_bond)
              bond_force.addBond(new_bond[0],new_bond[1],bond_length,cgmodel.bond_force_constant)
              if cgmodel.constrain_bonds:
               system.addConstraint(new_bond[0],new_bond[1], bond_length)
         system.addForce(bond_force)

        if cgmodel.include_nonbonded_forces:
         # Create nonbonded forces
         nonbonded_force = mm.NonbondedForce()
         bead_index = 0
         for monomer in range(cgmodel.polymer_length):
          for backbone_bead in range(cgmodel.backbone_length):
            system.addParticle(cgmodel.mass)
            nonbonded_force.addParticle(charge,sigma,epsilon)
            if backbone_bead in cgmodel.sidechain_positions:
              for sidechain_bead in range(cgmodel.sidechain_length):
                system.addParticle(cgmodel.mass)
                nonbonded_force.addParticle(charge,sigma,epsilon)
         system.addForce(nonbonded_force)
         nonbonded_force.createExceptionsFromBonds(new_bond_list,1.0,1.0)

        if cgmodel.include_bond_angle_forces:
         # Create bond angle potentials
         angle_list = cgmodel.get_angles()
         angle_force = mm.HarmonicAngleForce()
         for angle in angle_list:
              angle_force.addAngle(angle[0],angle[1],angle[2],cgmodel.equil_bond_angle,cgmodel.bond_angle_force_constant)
         system.addForce(angle_force)

        if cgmodel.include_torsion_forces:
         # Create torsion potentials
         torsion_list = cgmodel.get_dihedrals()
         torsion_force = mm.PeriodicTorsionForce()
         for torsion in torsion_list:
              torsion_force.addTorsion(torsion[0],torsion[1],torsion[2],torsion[3],1,cgmodel.equil_dihedral_angle,cgmodel.torsion_force_constant)
         system.addForce(torsion_force)
        
        return(system)

class CGModel(object):
        """
        Construct a coarse grained model.

        Parameters
        ----------

        positions: Positions for all of the particles, default = None

        polymer_length: Number of monomer units (integer), default = 8
      
        backbone_lengths: List of integers defining the umber of beads in the backbone for each monomer type
        portion of each (individual) monomer (integer), default = [1]

        sidechain_lengths: List of integers defining the umber of beads in the sidechain for each monomer type
        portion of each (individual) monomer (integer), default = [1]

        sidechain_positions: List of integers defining the backbone
        bead indices upon which we will place the sidechains,
        default = [0] (Place a sidechain on the backbone bead with
        index "0" (first backbone bead) in each (individual) monomer

        masses: Masses of all particle types
        ( List ( [ [ Backbone masses ], [ Sidechain masses ] ] ) )
        default = [ [ 12.0 * unit.amu ], [ 12.0 * unit.amu ] ]

        sigmas: Non-bonded bead Lennard-Jones equilibrium interaction distance
        ( [ [ float * simtk.unit.distance ], [ float * simtk.unit.distance ], [ float * simtk.unit.distance ] ] )
        default = [[8.4 * unit.angstrom],[8.4 * unit.angstrom],[8.4 * unit.angstrom]]

        epsilons: Non-bonded Lennard-Jones equilibrium interaction strengths
        ( [ [ float * simtk.unit.energy ], [ float * simtk.unit.energy ], [ float * simtk.unit.energy ] ] )
        default = [[0.5 * unit.kilocalorie_per_mole],[0.5 * unit.kilocalorie_per_mole],[0.5 * unit.kilocalorie_per_mole]]

        bond_lengths: Bond lengths for all bond types
        ( float * simtk.unit.distance )
        default = [[1.0 * unit.angstrom],[1.0 * unit.angstrom],[1.0 * unit.angstrom]]

        bond_force_constants: Bond force constants for all bond types
        ( float )
        default = [[9.9e5 kJ/mol/A^2],[9.9e5 kJ/mol/A^2],[9.9e5 kJ/mol/A^2]]

        charges: Charges for all beads
        ( float * simtk.unit.charge )
        default = [[0.0 * unit.elementary_charge],[0.0 * unit.elementary_charge]]

        num_beads: Total number of particles in the coarse grained model
        ( integer )
        default = polymer_length * ( backbone_length + sidechain_length )

        system: OpenMM system object, which stores forces, and can be used
        to check a model for energy conservation
        ( OpenMM System() class object )
        default = None

        topology: OpenMM topology object, which stores bonds, angles, and
        other structural attributes of the coarse grained model
        ( OpenMM Topology() class object )
        default = None

        constrain_bonds: Logical variable determining whether bond constraints
        are applied during a molecular dynamics simulation of the system.
        ( Logical )
        default = False

        include_bond_forces: Include contributions from bond
        (harmonic) potentials when calculating the potential energy
        ( Logical )
        default = True

        include_nonbonded_forces: Include contributions from nonbonded
        interactions when calculating the potential energy
        ( Logical )
        default = True

        include_bond_angle_forces: Include contributions from bond angles
        when calculating the potential energy
        ( Logical )
        default = False

        include_torsion_forces: Include contributions from torsions
        when calculating the potential energy
        ( Logical )
        default = False

        Attributes
        ----------

        polymer_length
        backbone_lengths
        sidechain_lengths
        sidechain_positions
        masses
        sigmas
        epsilons
        bond_lengths
        nonbonded_interaction_list
        bond_list
        bond_angle_list
        torsion_list
        bond_force_constants
        bond_angle_force_constants
        torsion_force_constants
        charges
        num_beads
        positions
        system
        topology
        constrain_bonds
        include_bond_forces
        include_nonbonded_forces
        include_bond_angle_forces
        include_torsion_forces

        Notes
        -----
        
        """

        # Built in class attributes
        _BUILT_IN_REGIONS = ('polymer_length','backbone_lengths','sidechain_lengths','sidechain_positions','masses','sigmas','epsilons','bond_lengths','bond_force_constants','bond_angle_force_constants','torsion_force_constants','charges','num_beads','positions','system','topology','constrain_bonds','bond_list','nonbonded_interaction_list','bond_angle_list','torsion_list','include_bond_forces','include_nonbonded_forces','include_bond_angle_forces','include_torsion_forces')

        def __init__(self, positions = None, polymer_length = 12, backbone_lengths = [1], sidechain_lengths = [1], sidechain_positions = [0], masses = [[12.0 * unit.amu],[12.0 * unit.amu]], sigmas = [[8.4 * unit.angstrom],[8.4 * unit.angstrom],[8.4 * unit.angstrom]], epsilons = [[0.5 * unit.kilocalorie_per_mole],[0.5 * unit.kilocalorie_per_mole],[0.5 * unit.kilocalorie_per_mole]], bond_lengths = [[1.0 * unit.angstrom],[1.0 * unit.angstrom],[1.0 * unit.angstrom]], bond_force_constants = [[9.9e5],[9.9e5],[9.9e5]], bond_angle_force_constants=[[200],[200],[200]],torsion_force_constants=[[200],[200],[200]], charges = [[0.0 * unit.elementary_charge],[0.0 * unit.elementary_charge]], constrain_bonds = False,include_bond_forces=True,include_nonbonded_forces=True,include_bond_angle_forces=True,include_torsion_forces=True,check_energy_conservation=True):

          """
          Initialize variables that were passed as input
          """

          self.polymer_length = polymer_length
          self.backbone_lengths = backbone_lengths
          self.sidechain_lengths = sidechain_lengths
          self.sidechain_positions = sidechain_positions
          self.monomer_types = self.get_monomer_types()
          self.num_beads = self.get_num_beads()
          self.masses = masses
          self.sigmas = sigmas
          self.epsilons = epsilons
          self.bond_lengths = bond_lengths
          self.bond_force_constants = bond_force_constants
          self.bond_angle_force_constants = bond_angle_force_constants
          self.torsion_force_constants = torsion_force_constants
          self.charges = charges

          if len(sigmas) != 0: include_nonbonded_forces = True
          if len(bond_force_constants) != 0: include_bond_forces = True
          if len(bond_angle_force_constants) != 0: include_bond_angle_forces = True
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
          list_of_masses = get_particle_masses(self)

          """
          Initialize new (coarse grained) particle types:
          """
          add_new_elements(self,list_of_masses)

          self.system = build_system(self)

          if positions == None: self.positions = util.random_positions(self) 
          else: self.positions = positions

        def get_num_beads(self):
          """
          Calculate the number of beads in our coarse grained model(s)
          """
          num_beads = []
          for backbone_length,sidechain_length in zip(self.backbone_lengths,self.sidechain_lengths):
           num_beads.append(self.polymer_length * ( backbone_length + sidechain_length ) )
          return(num_beads)

        def get_bond_list(self):
          """
          Construct a bond list for the coarse grained model
          """
          bond_list = []
          bead_index = 1
          for monomer in range(self.polymer_length):
            for backbone_bead in range(1,self.backbone_length+1):

             parent_index = get_parent_bead(self,bead_index,backbone_bead,sidechain_bead=False)
             if bead_index != 1:
              if parent_index < 0:
               print("Error: identifying parent index incorrectly when assigning a bond.")
               print("The bead index is: "+str(bead_index))
               print("The parent index is: "+str(parent_index))
               print("The backbone index is: "+str(backbone_bead))
               exit() 
              if parent_index != -1:
               if parent_index < bead_index:
                bond_list.append([parent_index,bead_index])
               else:
                bond_list.append([bead_index,parent_index])
             bead_index = bead_index + 1
             
             if backbone_bead-1 in self.sidechain_positions:
                for sidechain_bead in range(self.sidechain_length):
                  parent_index = get_parent_bead(self,bead_index,backbone_bead,sidechain_bead=True)
                  if parent_index < bead_index:
                   bond_list.append([parent_index,bead_index])
                  else:
                   bond_list.append([bead_index,parent_index])
                  bead_index = bead_index + 1

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
          Construct a list of bond angles for our coarse grained model
          """

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

          torsions = []
          for bond_1 in bond_list:
            torsion = [bond_1[0],bond_1[1]]
            for bond_2 in bond_list:
             if bond_2 != bond_1 and [bond_2[1],bond_2[0]] != bond_1:
              if bond_1[0] in bond_2 or bond_1[1] in bond_2:
               if bond_2[0] not in bond_angle:
                torsion.append(bond_2[0])
               if bond_2[1] not in bond_angle:
                bond_angle.append(bond_2[1])

             for bond_3 in bond_list:
              if bond_3 != bond_1 and [bond_3[1],bond_3[0]] != bond_1:
                if bond_3 != bond_2 and [bond_3[1],bond_3[0]] != bond_2:
                  if bond_3[0] in bond_2 or bond_1[1] in bond_2:
                    if bond_2[0] not in bond_angle:
                      torsion.append(bond_2[0])
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

          return(dihedrals)
