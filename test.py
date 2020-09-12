import numpy as np  
import pandas as pd

flag_shape = 0 
flag_count = 0

with open("2gktH.pdb", "r") as pdb_file: 
    for line in pdb_file:
        if line.startswith("ATOM"): 
            flag_shape += 1
               
p = np.full((flag_shape,9),np.nan) 
table = pd.DataFrame(p)
table.columns = ["AA","position_AA","atom","position_atom","chain","coord_x","coord_y","coord_z","occupancy"]    

with open("2gktH.pdb", "r") as pdb_file: 
    for line in pdb_file:
        if line.startswith("ATOM"): 
            AA_name = line[17:20].strip() 
            AA_position = line[22:26].strip()
            atom_name = line[12:16].strip() 
            atom_position = line[6:11].strip()
            chain = line[21]
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54]) 
            occupancy = line[54:60].strip()   
            table.iloc[flag_count] = [AA_name,AA_position,atom_name,atom_position,chain,x,y,z,occupancy]
            flag_count += 1
            
table['position_atom'] = table['position_atom'].replace(['0'],np.nan)
print(table)
