def lecture(name):
    flag_count = 0
    flag_shape = 0 
    with open(name, "r") as pdb_file: 
        for line in pdb_file:
            if line.startswith("ATOM"): 
                flag_shape += 1   
    p = np.full((flag_shape,9),np.nan) 
    table = pd.DataFrame(p)
    table.columns = ["AA","position_AA","atom","position_atom","chain","coord_x","coord_y","coord_z","occupancy"]    
    with open(name, "r") as pdb_file: 
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
        return(table) 

def creation_table(table): 
    table_sulfuric = pd.DataFrame((table.loc[(table['AA'] == "CYS") | (table['AA'] == "MET")])) 
    table_hydrophobic = pd.DataFrame((table.loc[(table['AA'] == "ALA") | (table['AA'] == "VAL") | \
    					(table['AA'] == "LEU") | (table['AA'] == "ILE") | (table['AA'] == "MET") | \
    					(table['AA'] == "PHE") | (table['AA'] == "TRP") | (table['AA'] == "PRO") | \
    					(table['AA'] == "TYR")]))
    return(table_sulfuric, table_hydrophobic)

def calcule_distance(ss_table):
	i = 0
	j = 0
	liste_distance = []
	for i in range(len(ss_table) -1):
		j = j+i
		for j in range(i+1, len(ss_table)):
			distance = (np.sqrt(((ss_table[i,5] - ss_table[(j),5])**2) + \
						((ss_table[i,6] - ss_table[(j),6])**2) + \
						((ss_table[i,7] - ss_table[(j),7])**2)))
			liste_distance.append(distance) 
	
	return(liste_distance)



if __name__ == "__main__": 
    import numpy as np
    import pandas as pd 
    var_lecture = lecture("2gktH.pdb") 
    var_creation_table = creation_table(var_lecture)
    ss_table_sulfuric = np.array(var_creation_table[0].loc[(var_creation_table[0]["atom"] == "SG")])
    distance_sulfuric = calcule_distance(ss_table_sulfuric)
    print(distance_sulfuric)