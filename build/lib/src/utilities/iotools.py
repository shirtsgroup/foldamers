# This script generates random coordinates for a CG polymer

# =============================================================================================
# 1) PYTHON PACKAGE IMPORTS
# =============================================================================================

# System packages
from simtk.openmm.app.pdbfile import PDBFile
import numpy as np
import math, random
from simtk import unit
# =============================================================================================
# 2) ENVIRONMENT/JOB SETTINGS
# =============================================================================================

def write_bonds(CGModel,pdb_object):
        bead_index = 1
        for monomer_index in range(CGModel.polymer_length):
          for backbone_bead in range(CGModel.backbone_length):

            if backbone_bead != 0:
             if backbone_bead - 1 in CGModel.sidechain_positions:
              parent_bead = bead_index - CGModel.sidechain_length - 1
             else:
              parent_bead = bead_index - 1

             pdb_object.write("CONECT"+str("{:>5}".format(bead_index))+str("{:>5}".format(parent_bead))+"\n")      
             bead_index = bead_index + 1

            else:
             if bead_index != 1:
              if backbone_bead - 1 in CGModel.sidechain_positions:
               parent_bead = bead_index - CGModel.sidechain_length - 1
              else:
               parent_bead = bead_index - 1
              pdb_object.write("CONECT"+str("{:>5}".format(bead_index))+str("{:>5}".format(parent_bead))+"\n")
              bead_index = bead_index + 1

             if bead_index == 1:
              bead_index = bead_index + 1
            if backbone_bead in CGModel.sidechain_positions:
              for sidechain_bead in range(0,CGModel.sidechain_length):
               pdb_object.write("CONECT"+str("{:>5}".format(bead_index))+str("{:>5}".format(bead_index-1))+"\n")
               bead_index = bead_index + 1
        pdb_object.write(str("END\n"))
        return

def write_cg_pdb(cgmodel,topology,positions,file_name):
        file_obj = open(file_name,'w')
        PDBFile.writeHeader(topology, file_obj)
        PDBFile.writeModel(topology, positions, file_obj)
        write_bonds(cgmodel,file_obj)
        file_obj.close()
        return


def write_pdbfile(CGModel,filename):
        """
        Writes the positions in 'CGModel' to the file 'filename'.

        Parameters
        ----------

        CGModel: Coarse grained model class object

        filename: Path to the file where we will write PDB coordinates.

        """

        pdb_object = open(filename,"w")


        coordinates = CGModel.positions
        bead_index = 1
        for monomer_index in range(CGModel.polymer_length):
          element_index = 1
          for backbone_bead in range(CGModel.backbone_length):

            if monomer_index in list([0,CGModel.polymer_length-1]):
             pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"   B"+str(element_index)+str("{:>4}".format(str("MT")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
            else:
             pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"   B"+str(element_index)+str("{:>4}".format(str("M")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
            bead_index = bead_index + 1
            element_index = element_index + 1

            if backbone_bead in CGModel.sidechain_positions:
              for sidechain_bead in range(0,CGModel.sidechain_length):
                if monomer_index in list([0,CGModel.polymer_length-1]):
                 pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"   S"+str(element_index)+str("{:>4}".format(str("MT")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
                else:
                 pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"   S"+str(element_index)+str("{:>4}".format(str("M")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
                bead_index = bead_index + 1
                element_index = element_index + 1
        pdb_object.write(str("TER\n"))

        write_bonds(CGModel,pdb_object)
        pdb_object.close()
        return
