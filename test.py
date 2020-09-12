import numpy as np 
table = np.array([["AA","Position of AA","Atom","Position of Atom","Chain","Xcoord","Ycoord","Zcoord","Occupancy"]])
with open("2gkt.pdb", "r") as pdb_file:
    for line in pdb_file:
        if line.startswith("ATOM"): 
            AA_name = line[17:20].strip() 
            AA_position = line[22:26].strip()
            atome_name = line[12:16].strip() 
            atome_position = line[6:11].strip()
            chaine = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54]) 
            occupancy = line[54:60].strip()  
print(table)