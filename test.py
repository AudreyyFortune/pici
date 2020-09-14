""" 
Projet Court.
Sujet 6 : Protein Interactions Calculator implémentation.

----------------------------------------------------------------
Nom_1 : BATISTA Reda
Nom_2 : FORTUNE Audrey
Cursus : M2 Bio informatique
----------------------------------------------------------------
"""



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
    df_disulfure = pd.DataFrame((table.loc[(table['AA'] == "CYS") & (table['atom'] == "SG")]))
    disulfure = np.array(df_disulfure)
    #table_hydrophobic = pd.DataFrame((table.loc[(table['AA'] == "ALA") | (table['AA'] == "VAL") | \
    					#(table['AA'] == "LEU") | (table['AA'] == "ILE") | (table['AA'] == "MET") | \
    					#(table['AA'] == "PHE") | (table['AA'] == "TRP") | (table['AA'] == "PRO") | \
    					#(table['AA'] == "TYR")]))
    df_ionic_coo = pd.DataFrame((table.loc[(((table['AA'] == "ASP") & (table['atom'] == "OD1")) | \
    			(table['AA'] == "GLU") & (table['atom'] == "OE1"))]))
    df_ionic_nh = pd.DataFrame((table.loc[((((table['AA'] == "ARG") & (table['atom'] == "NH1")) | \
    			(table['AA'] == "LYS") & (table['atom'] == "NZ")) | ((table['AA'] == "HIS") & (table['atom'] == "ND1")))]))
    return(disulfure, df_ionic_coo, df_ionic_nh)#, table_hydrophobic)

def calcule_distance(ss_table):
	i = 0
	j = 0
	compteur = 0
	liste_distance = []
	for i in range(len(ss_table) -1):
		j = j+i
		for j in range(i+1, len(ss_table)):
			distance = (np.sqrt(((ss_table[i,5] - ss_table[(j),5])**2) + ((ss_table[i,6] \
						- ss_table[(j),6])**2) + ((ss_table[i,7] - ss_table[(j),7])**2)))
			liste_distance.append(ss_table[i,1])
			liste_distance.append(ss_table[i,0])
			liste_distance.append(ss_table[i,4])
			liste_distance.append(ss_table[j,1])
			liste_distance.append(ss_table[j,0])
			liste_distance.append(ss_table[j,4])			
			liste_distance.append(distance) 
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue', \
				 'chain', 'distance'], data = array)
	array['distance'] = array['distance'].astype(float)
	return(array)

def maybe(t1, t2):
	i = 0
	j = 0
	compteur = 0
	liste_distance = []
	for i in range(len(t1)):
		j += 1
		for j in range(i, len(t2)):
			distance = (np.sqrt(((t1[i,5] - t2[(j),5])**2) + ((t1[i,6] - t2[(j),6])**2) + ((t1[i,7] - t2[(j),7])**2)))
			liste_distance.append(t1[i,1])
			liste_distance.append(t1[i,0])
			liste_distance.append(t1[i,4])
			liste_distance.append(t2[j,1])
			liste_distance.append(t2[j,0])
			liste_distance.append(t2[j,4])			
			liste_distance.append(distance) 
			compteur += 1
		return(liste_distance)
		


def cutoff(array):
    table_disulfure = array[array["distance"] < 2.2]	
    return(table_disulfure)


if __name__ == "__main__": 
    import sys
    import os
    import numpy as np
    import pandas as pd
    if len(sys.argv) != 1:
    	sys.exit("Erreur : Il faut exactement un argument.")
    Question = input("Veuillez entrer le nom du fichier :")
    Reponse = str(Question)
    comp_line = 0
    if os.path.exists(Reponse):
    	print("Le fichier est présent.")
    else:
    	sys.exit("Le fichier est absent.")
    with open(Reponse,"r") as pdb_files:
    	for lines in pdb_files:
    		if lines.startswith("ATOM"):
    			comp_line=comp_line+1
    	if comp_line > 3:
    		var_lecture = lecture(Reponse)
    	#var_lecture = lecture("2gktH.pdb")
    		var_creation_table = creation_table(var_lecture)
    		distance_disulfure = calcule_distance(var_creation_table[0])
    		print(var_creation_table)
    		distance_ionic = maybe(var_creation_table[1],var_creation_table[2])
    		print(distance_ionic)
    		table_disulfure = cutoff(distance_disulfure)
    		print("""Intraprotein Disulphide Bridges \n Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms \n""", table_disulfure)
