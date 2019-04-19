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
          for backbone_bead in range(CGModel.backbone_length):

            if monomer_index in list([0,CGModel.polymer_length-1]):
             pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  CG1"+str("{:>4}".format(str("MT")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
            else:
             pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  CG1"+str("{:>4}".format(str("M")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
            bead_index = bead_index + 1

            if backbone_bead in CGModel.sidechain_positions:
              for sidechain_bead in range(0,CGModel.sidechain_length):
                if monomer_index in list([0,CGModel.polymer_length-1]):
                 pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  CG2"+str("{:>4}".format(str("MT")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
                else:
                 pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  CG2"+str("{:>4}".format(str("M")))+" A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1]._value,3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2]._value,3)))+"  1.00  0.00\n"))
                bead_index = bead_index + 1
        pdb_object.write(str("TER\n"))

        bead_index = 1
        for monomer_index in range(CGModel.polymer_length):
          for backbone_bead in range(CGModel.backbone_length):

            if bead_index != 1:
             parent_bead = bead_index - CGModel.sidechain_length - 1
             pdb_object.write("CONECT"+str("{:>5}".format(bead_index))+str("{:>5}".format(parent_bead))+"\n")            
             bead_index = bead_index + 1

            if bead_index == 1:
             bead_index = bead_index + 1

            if backbone_bead in CGModel.sidechain_positions:
              for sidechain_bead in range(0,CGModel.sidechain_length):
                pdb_object.write("CONECT"+str("{:>5}".format(bead_index))+str("{:>5}".format(bead_index-1))+"\n")
                bead_index = bead_index + 1
        pdb_object.write(str("END\n"))
        pdb_object.close()
        return
