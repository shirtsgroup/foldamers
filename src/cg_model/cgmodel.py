from simtk import unit
import sys, os
from foldamers.src.utilities import util
from simtk import openmm as mm
from simtk.openmm.app.topology import Topology
from simtk.openmm.app.topology import Residue
import simtk.openmm.app.element as elem
from itertools import chain, combinations, product

def get_particle_type(cgmodel,particle_index,particle_name=None):
        """
        Returns the name of a particle, given its index within the model

        Parameters
        ----------

        cgmodel: CGModel() class object

        particle_index: Index of the particle for which we would like to determine the type
        Type: int()

        Returns
        -------

        particle_type: 'backbone' or 'sidechain'
        Type: str()

        """
        if particle_name == None: particle_name = cgmodel.particle_list[particle_index]
        if 'B' in particle_name: particle_type = 'backbone'
        if 'S' in particle_name: particle_type = 'sidechain'

        return(particle_type)

def get_particle_mass(cgmodel,particle_index):
        """
        Returns the mass for a particle, given its index.

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        Mass

        """
        particle_type = get_particle_type(cgmodel,particle_index)
        if particle_type == 'backbone': particle_mass = cgmodel.masses['backbone_bead_masses']
        if particle_type == 'sidechain': particle_mass = cgmodel.masses['sidechain_bead_masses']
        return(particle_mass)

def get_particle_charge(cgmodel,particle_index):
        """
        Returns the charge for a particle, given its index.

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        Charge

        """
        particle_type = get_particle_type(cgmodel,particle_index)
        if particle_type == 'backbone': particle_charge = cgmodel.charges['backbone_bead_charges']
        if particle_type == 'sidechain': particle_charge = cgmodel.charges['sidechain_bead_charges']
        return(particle_charge)

def get_sigma(cgmodel,particle_index,particle_type=None):
        """
        Returns the sigma value for a particle, given its index within the coarse grained model.

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        Sigma

        """

        if particle_type == None: particle_type = get_particle_type(cgmodel,particle_index)
        if particle_type == 'backbone': sigma = cgmodel.sigmas['bb_bb_sigma']
        if particle_type == 'sidechain': sigma = cgmodel.sigmas['sc_sc_sigma']
        return(sigma)

def get_epsilon(cgmodel,particle_index,particle_type=None):
        """
        Returns the epsilon value for a particle, given its index.

        Parameters
        ----------

        cgmodel: CGModel() class object

        Returns
        -------

        Epsilon

        """
        if particle_type == None: particle_type = get_particle_type(cgmodel,particle_index)
        if particle_type == 'backbone': epsilon = cgmodel.epsilons['bb_bb_eps']
        if particle_type == 'sidechain': epsilon = cgmodel.epsilons['sc_sc_eps']
        return(epsilon)

def get_all_particle_masses(cgmodel):
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
        list_of_masses.append(cgmodel.masses['backbone_bead_masses'])
        list_of_masses.append(cgmodel.masses['sidechain_bead_masses'])
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
        particle_list = []
        for monomer_type in cgmodel.monomer_types:
         for backbone_bead in range(monomer_type['backbone_length']):
          particle_name = str("bb-"+str(cg_particle_index))
          particle_symbol = str("B"+str(cg_particle_index))
          if particle_symbol not in elem.Element._elements_by_symbol:
           elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
           particle_list.append(particle_symbol)
           element_index = element_index + 1
           cg_particle_index = cg_particle_index + 1
           mass_index = mass_index + 1
          if type(monomer_type['sidechain_positions']) == int:
           sidechain_positions = [monomer_type['sidechain_positions']]
          else:
           sidechain_positions = monomer_type['sidechain_positions']
          if backbone_bead in sidechain_positions:
           for sidechain in range(monomer_type['sidechain_length']):
            if particle_symbol not in elem.Element._elements_by_symbol:
             particle_name = str("sc-"+str(cg_particle_index))
             particle_symbol = str("S"+str(cg_particle_index))
             elem.Element(element_index,particle_name,particle_symbol,list_of_masses[mass_index])
             particle_list.append(particle_symbol)
             element_index = element_index + 1
             cg_particle_index = cg_particle_index + 1
             mass_index = mass_index + 1
        return(particle_list)

def get_bond_length_from_names(cgmodel,particle_1_name,particle_2_name):
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

        if 'B' in particle_2_name: particle_2_type = 'backbone'
        else: particle_2_type = 'sidechain'

        if particle_1_type == 'backbone' and particle_2_type == 'backbone':
         bond_length = cgmodel.bond_lengths['bb_bb_bond_length']

        if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
         bond_length = cgmodel.bond_lengths['bb_sc_bond_length']

        if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
         bond_length = cgmodel.bond_lengths['bb_bb_bond_length']

        return(bond_length)


def get_bond_length(cgmodel,particle_1_index,particle_2_index):
        """
        Determines the correct bond force constant for two particles

        Parameters
        ----------

        cgmodel: CGModel() class object

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
        if 'B' in cgmodel.particle_list[particle_1_index]: particle_1_type = 'backbone'
        else: particle_1_type = 'sidechain'

        if 'B' in cgmodel.particle_list[particle_2_index]: particle_2_type = 'backbone'
        else: particle_2_type = 'sidechain'

        if particle_1_type == 'backbone' and particle_2_type == 'backbone':
         bond_length = cgmodel.bond_lengths['bb_bb_bond_length']

        if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
         bond_length = cgmodel.bond_lengths['bb_sc_bond_length']

        if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
         bond_length = cgmodel.bond_lengths['bb_bb_bond_length']

        return(bond_length)

def get_bond_force_constant(cgmodel,particle_1_index,particle_2_index):
        """
        Determines the correct bond force constant for two particles

        Parameters
        ----------

        cgmodel: CGModel() class object

        particle_1_index: Index of the first particle in the bond
        ( integer )
        Default = None

        particle_2_index: Index of the second particle in the bond
        ( integer )
        Default = None

        Returns
        -------

        bond_force_constant: Bond force constant for the bond defined by these two particles.
        ( Integer )

        """
        if 'B' in cgmodel.particle_list[particle_1_index]: particle_1_type = 'backbone'
        else: particle_1_type = 'sidechain'

        if 'B' in cgmodel.particle_list[particle_2_index]: particle_2_type = 'backbone'
        else: particle_2_type = 'sidechain'

        if particle_1_type == 'backbone' and particle_2_type == 'backbone':
         bond_force_constant = cgmodel.bond_force_constants['bb_bb_bond_k']

        if particle_1_type == 'backbone' and particle_2_type == 'sidechain':
         bond_force_constant = cgmodel.bond_force_constants['bb_sc_bond_k']

        if particle_1_type == 'sidechain' and particle_2_type == 'sidechain':
         bond_force_constant = cgmodel.bond_force_constants['bb_bb_bond_k']

        return(bond_force_constant)

def get_torsion_force_constant(cgmodel,torsion):
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
        if 'B' in cgmodel.particle_list[torsion[0]]: particle_types[0] = 'backbone'
        else: particle_types[0] = 'sidechain'

        if 'B' in cgmodel.particle_list[torsion[1]]: particle_types[1] = 'backbone'
        else: particle_types[1] = 'sidechain'

        if 'B' in cgmodel.particle_list[torsion[2]]: particle_types[2] = 'backbone'
        else: particle_types[2] = 'sidechain'

        if 'B' in cgmodel.particle_list[torsion[3]]: particle_types[3] = 'backbone'
        else: particle_types[3] = 'sidechain'

        if particle_types[0] == 'sidechain':
         if particle_types[1] == 'backbone':
          if particle_types[2] == 'backbone':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_bb_bb_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_bb_bb_sc_torsion_k']
          if particle_types[2] == 'sidechain':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_sc_sc_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_bb_sc_sc_torsion_k']
         if particle_types[1] == 'sidechain':
          if particle_types[2] == 'backbone':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_sc_bb_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_sc_bb_sc_torsion_k']
          if particle_types[2] == 'sidechain':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_sc_sc_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['sc_sc_sc_sc_torsion_k']
        if particle_types[0] == 'backbone':
         if particle_types[1] == 'backbone':
          if particle_types[2] == 'backbone':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_bb_bb_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_bb_bb_sc_torsion_k']
          if particle_types[2] == 'sidechain':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_bb_sc_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_bb_sc_sc_torsion_k']
         if particle_types[1] == 'sidechain':
          if particle_types[2] == 'backbone':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_sc_bb_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_sc_bb_sc_torsion_k']
          if particle_types[2] == 'sidechain':
           if particle_types[3] == 'backbone':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_sc_sc_bb_torsion_k']
           if particle_types[3] == 'sidechain':
            torsion_force_constant = cgmodel.torsion_force_constants['bb_sc_sc_sc_torsion_k']
        return(torsion_force_constant)

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

def build_topology(cgmodel):
        """

        Construct an OpenMM topology for our coarse grained model

        Parameters
        ----------

        polymer_length: Number of monomers in our coarse grained model
        ( integer )

        backbone_length: Number of backbone beads on individual monomers
        in our coarse grained model, ( integer )

        sidechain_length: Number of sidechain beads on individual monomers
        in our coarse grained model, ( integer )

        """

        topology = Topology()

        chain = topology.addChain()
        residue_index = 1
        cg_particle_index = 1
        for monomer_type in cgmodel.sequence:
         residue = topology.addResidue(str(residue_index), chain)
         for backbone_bead in range(monomer_type['backbone_length']): 
          particle_name = str("bb-"+str(cg_particle_index))
          particle_symbol = str("B"+str(cg_particle_index))
          particle = topology.addAtom(particle_symbol, particle_name, residue)
          if backbone_bead == 0 and residue_index != 1:
           topology.addBond(particle,last_backbone_particle)
          last_backbone_particle = particle
          cg_particle_index = cg_particle_index + 1
          if backbone_bead in [monomer_type['sidechain_positions']]:
           for sidechain_bead in range(monomer_type['sidechain_length']):
             particle_name = str("sc-"+str(cg_particle_index))
             particle_symbol = str("S"+str(cg_particle_index))
             particle = topology.addAtom(particle_symbol, particle_name, residue)
             if sidechain_bead == 0:
              topology.addBond(particle,last_backbone_particle)
             if sidechain_bead != 0:
              topology.addBond(particle,last_sidechain_particle)
             last_sidechain_particle = particle
             cg_particle_index = cg_particle_index + 1
         residue_index = residue_index + 1
        return(topology)



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
#        sigma = cgmodel.sigma.in_units_of(unit.nanometer)._value
#        charge = cgmodel.charge._value
#        epsilon = cgmodel.epsilon.in_units_of(unit.kilojoule_per_mole)._value
#        bond_length = cgmodel.bond_length.in_units_of(unit.nanometer)._value

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
              bond_force_constant = get_bond_force_constant(cgmodel,new_bond[0],new_bond[1])
              bond_length = get_bond_length(cgmodel,new_bond[0],new_bond[1])
              bond_length = bond_length.in_units_of(unit.nanometer)._value
              bond_force.addBond(new_bond[0],new_bond[1],bond_length,bond_force_constant)
              if cgmodel.constrain_bonds:
               system.addConstraint(new_bond[0],new_bond[1], bond_length)
         system.addForce(bond_force)

        if cgmodel.include_nonbonded_forces:
         # Create nonbonded forces
         nonbonded_force = mm.NonbondedForce()
         bead_index = 0
         for monomer_type in cgmodel.sequence:
          for backbone_bead in range(monomer_type['backbone_length']):
            mass = get_particle_mass(cgmodel,bead_index)
            charge = get_particle_charge(cgmodel,bead_index)
            sigma = get_sigma(cgmodel,bead_index)
            epsilon = get_epsilon(cgmodel,bead_index)
            system.addParticle(mass)
            bead_index = bead_index + 1
            sigma = sigma.in_units_of(unit.nanometer)._value
            charge = charge._value
            epsilon = epsilon.in_units_of(unit.kilojoule_per_mole)._value
            nonbonded_force.addParticle(charge,sigma,epsilon)
            if backbone_bead in [monomer_type['sidechain_positions']]:
              for sidechain_bead in range(monomer_type['sidechain_length']):
                mass = get_particle_mass(cgmodel,bead_index)
                charge = get_particle_charge(cgmodel,bead_index)
                sigma = get_sigma(cgmodel,bead_index)
                epsilon = get_epsilon(cgmodel,bead_index)
                system.addParticle(mass)
                bead_index = bead_index + 1
                sigma = sigma.in_units_of(unit.nanometer)._value
                charge = charge._value
                epsilon = epsilon.in_units_of(unit.kilojoule_per_mole)._value
                nonbonded_force.addParticle(charge,sigma,epsilon)
         system.addForce(nonbonded_force)
         nonbonded_force.createExceptionsFromBonds(new_bond_list,1.0,1.0)

        if cgmodel.include_bond_angle_forces:
         # Create bond angle potentials
         angle_list = cgmodel.get_bond_angle_list()
         angle_force = mm.HarmonicAngleForce()
         for angle in angle_list:
              angle_force.addAngle(angle[0],angle[1],angle[2],cgmodel.equil_bond_angle,cgmodel.bond_angle_force_constant)
         system.addForce(angle_force)

        if cgmodel.include_torsion_forces:
         # Create torsion potentials
         torsion_list = cgmodel.get_torsion_list()
         torsion_force = mm.PeriodicTorsionForce()
         for torsion in torsion_list:
              torsion_force_constant = get_torsion_force_constant(cgmodel,torsion)
              torsion_force.addTorsion(torsion[0],torsion[1],torsion[2],torsion[3],1,cgmodel.equil_dihedral_angle,torsion_force_constant)
         system.addForce(torsion_force)
        
        return(system)


def basic_cgmodel(polymer_length=8,backbone_length=1,sidechain_length=1,sidechain_positions=[0],mass=12.0 * unit.amu,charge=0.0 * unit.elementary_charge,bond_length=1.0 * unit.angstrom,sigma=2.5*unit.angstrom,epsilon=0.5 * unit.kilocalorie_per_mole,positions=None):
        """
        Given a minimal set of model parameters, this function creates a cgmodel class object.

        Parameters
        ----------

        polymer_length: Number of monomer units (integer)
        default = 8

        backbone_length: Integer defining the number of beads in the backbone
        default = 1

        sidechain_length: Integer defining the number of beads in the sidechain
        default = 1

        polymer_length: Number of monomer units (integer), default = 8

        sidechain_positions: List of integers defining the backbone
        bead indices upon which we will place the sidechains,
        default = [0]

        mass: Mass for all coarse grained beads.
        default = 12.0 * unit.amu

        bond_length: Bond length for all bond types
        Default = 1.0 * unit.angstrom

        sigma: Non-bonded bead Lennard-Jones equilibrium interaction distance.
        default = 2.5 * bond_length

        epsilon: Non-bonded Lennard-Jones equilibrium interaction energy
        default = 0.5 * unit.kilocalorie_per_mole

        charge: Charge for all particles
        default = 0.0 * unit.elementary_charge

        positions: Positions for coarse grained particles in the model.
        default = None

        Returns
        -------

        cgmodel: CGModel() class object

        """
        bond_force_constant = 9.9e5
        torsion_force_constant = 200
        equil_dihedral_angle = 180
        bond_angle_force_constant = 200
        backbone_lengths = [1] # Number of backbone beads in unique monomer types
        sidechain_lengths = [1] # Number of sidechain beads in unique monomer types
        sidechain_positions = [0] # Index of the backbone bead(s) to which sidechains are bonded
        polymer_length = 8 # Number of monomers in the polymer
        masses = {'backbone_bead_masses': mass, 'sidechain_bead_masses': mass} # List of bead masses
        sigmas = {'bb_bb_sigma': sigma,'bb_sc_sigma': sigma,'sc_sc_sigma': sigma} # Lennard-Jones interaction distances.  List of unique interaction types
        bond_lengths = {'bb_bb_bond_length': bond_length,'bb_sc_bond_length': bond_length,'sc_sc_bond_length': bond_length} # bond length
        bond_force_constants = {'bb_bb_bond_k': bond_force_constant,'bb_sc_bond_k': bond_force_constant, 'sc_sc_bond_k': bond_force_constant} # Units = kJ/mol/A^2 List of bond force constants for unique bond types
        epsilons = {'bb_bb_eps': epsilon,'bb_sc_eps': epsilon,'sc_sc_eps': epsilon} # Lennard-Jones interaction strength.  List of unique interaction types
        charges = {'backbone_bead_charges': charge,'sidechain_bead_charges': charge} # Charge of beads.
        torsion_force_constants = {'bb_bb_bb_bb_torsion_k': torsion_force_constant,'bb_bb_bb_sc_torsion_k': torsion_force_constant,'bb_bb_sc_sc_torsion_k': torsion_force_constant, 'bb_sc_sc_sc_torsion_k': torsion_force_constant, 'sc_bb_bb_sc_torsion_k': torsion_force_constant, 'bb_sc_sc_bb_torsion_k': torsion_force_constant, 'sc_sc_sc_sc_torsion_k': torsion_force_constant} # List of torsion force constants (k) for each of the torsions in our coarse grained model
        bond_angle_force_constants = {'bb_bb_bb_angle_k': bond_angle_force_constant,'bb_bb_sc_angle_k': bond_angle_force_constant,'bb_sc_sc_angle_k': bond_angle_force_constant,'sc_sc_sc_angle_k': bond_angle_force_constant} # List of bond angle force constants (k) for each of the bond angles in our coarse grained model
        if len(sigmas) != 0: include_nonbonded_forces = True
        if len(bond_force_constants) != 0: include_bond_forces = True
        # include_bond_forces = False
        if len(bond_angle_force_constants) != 0: include_bond_angle_forces = True
        if len(torsion_force_constants) != 0: include_torsion_forces = True
        include_bond_angle_forces = False
        cgmodel = CGModel(positions=positions,polymer_length=polymer_length,backbone_lengths=backbone_lengths, sidechain_lengths=sidechain_lengths, sidechain_positions = sidechain_positions, masses = masses, sigmas = sigmas, epsilons = epsilons, bond_lengths = bond_lengths, bond_force_constants = bond_force_constants, torsion_force_constants=torsion_force_constants, equil_dihedral_angle=equil_dihedral_angle,bond_angle_force_constants=bond_angle_force_constants, charges =charges,include_bond_forces=include_bond_forces,include_nonbonded_forces=include_nonbonded_forces,include_bond_angle_forces=include_bond_angle_forces,include_torsion_forces=include_torsion_forces,check_energy_conservation=False)
        return(cgmodel)


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
        equil_dihedral_angle
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
        _BUILT_IN_REGIONS = ('polymer_length','backbone_lengths','sidechain_lengths','sidechain_positions','masses','sigmas','epsilons','bond_lengths','bond_force_constants','bond_angle_force_constants','torsion_force_constants','equil_dihedral_angle','charges','num_beads','positions','system','topology','constrain_bonds','bond_list','nonbonded_interaction_list','bond_angle_list','torsion_list','include_bond_forces','include_nonbonded_forces','include_bond_angle_forces','include_torsion_forces')

        def __init__(self, positions = None, polymer_length = 12, backbone_lengths = [1], sidechain_lengths = [1], sidechain_positions = [0], masses = {'backbone_bead_masses': 12.0 * unit.amu, 'sidechain_bead_masses': 12.0 * unit.amu}, sigmas = {'bb_bb_sigma': 8.4 * unit.angstrom,'bb_sc_sigma': 8.4 * unit.angstrom,'sc_sc_sigma': 8.4 * unit.angstrom}, epsilons = {'bb_bb_eps': 0.5 * unit.kilocalorie_per_mole,'bb_sc_eps': 0.5 * unit.kilocalorie_per_mole,'sc_sc_eps': 0.5 * unit.kilocalorie_per_mole}, bond_lengths = {'bb_bb_bond_length': 1.0 * unit.angstrom,'bb_sc_bond_length': 1.0 * unit.angstrom,'sc_sc_bond_length': 1.0 * unit.angstrom}, bond_force_constants = {'bb_bb_bond_k': 9.9e5,'bb_sc_bond_k': 9.9e5, 'sc_sc_bond_k': 9.9e5}, bond_angle_force_constants={'bb_bb_bb_angle_k': 200,'bb_bb_sc_angle_k': 200,'bb_sc_sc_angle_k': 200,'sc_sc_sc_angle_k': 200}, torsion_force_constants={'bb_bb_bb_bb_torsion_k': 200,'bb_bb_bb_sc_torsion_k': 200,'bb_bb_sc_sc_torsion_k': 200, 'bb_sc_sc_sc_torsion_k': 200, 'sc_bb_bb_sc_torsion_k': 200, 'bb_sc_sc_bb_torsion_k': 200, 'sc_sc_sc_sc_torsion_k': 200}, equil_dihedral_angle = 180, charges = {'backbone_bead_charges': 0.0 * unit.elementary_charge,'sidechain_bead_charges': 0.0 * unit.elementary_charge}, constrain_bonds = False,include_bond_forces=True,include_nonbonded_forces=True,include_bond_angle_forces=True,include_torsion_forces=True,check_energy_conservation=True,homopolymer=True):

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
          list_of_masses = get_all_particle_masses(self)

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
