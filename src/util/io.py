
def write_positions_to_pdbfile(coordinates,filename,model_settings):
 box_size,polymer_length,backbone_length,sidechain_length,sidechain_positions = model_settings[:]
 monomer_size = backbone_length + sidechain_length
 pdb_object = open(filename,"w")
 bead_index = 1
 coordinates = coordinates._value
 for monomer_index in range(0,polymer_length):
  for backbone_bead in range(0,backbone_length):
   pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  X   CG  A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0],3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1],3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2],3)))+"  1.00  0.00\n"))
   bead_index = bead_index + 1
  for sidechain_bead in range(0,sidechain_length):
   pdb_object.write(str("ATOM"+str("{:>7}".format(bead_index))+"  Q   CG  A"+str("{:>4}".format(monomer_index+1))+"     "+str("{:>7}".format(round(coordinates[bead_index-1][0],3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][1],3)))+" "+str("{:>7}".format(round(coordinates[bead_index-1][2],3)))+"  1.00  0.00\n"))
   bead_index = bead_index + 1
 pdb_object.write(str("END"))
 pdb_object.close()
 return

def write_positions_to_xyzfile(coordinates,filename,model_settings):
 box_size,polymer_length,backbone_length,sidechain_length,sidechain_positions = model_settings[:]
 monomer_size = backbone_length + sidechain_length
 xyz_object = open(filename,"w")
 xyz_object.write(str(polymer_length * monomer_size )+"\n")
 xyz_object.write("\n")
 polymer_index = 1
 bead_index = 1
 while polymer_index <= polymer_length:
  monomer_index = 1
  while monomer_index <= monomer_size:
   xyz_object.write(str("C "+str("{:10.5f}".format(coordinates[bead_index-1][0]))+" "+str("{:10.5f}".format(coordinates[bead_index-1][1]))+" "+str("{:10.5f}".format(coordinates[bead_index-1][2]))+"\n")
)
   bead_index = bead_index + 1
   monomer_index = monomer_index + 1
#   if bead_index <= 9:
  polymer_index = polymer_index + 1
 xyz_object.close()
 return
