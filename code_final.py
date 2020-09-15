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
	"""
	Lecture d'un fichier pdb et extraction des données.

		Parameters 
		---------- 
		name : fichier pdb 
    		Contenant les coordonées atomiques des atomes des acides aminés.

		Returns
		------- 
		dataframe 
    		Contenant toutes les lignes du fichier pdb contenant les AA.
    		présenté en un tableau avec comme colonnes : le résidu, la position, l'atome,
    			la position de l'atome, la chaine, les coordonnées x, y et z, ainsi que l'occupancy.
	"""
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


def disulphide_bridges(table):
	"""
	Calcul des distances entre chaque atome SG des cystéines, et sélection des ponts disulfures.
		
		Parameters
		----------
			table : dataframe obtenu avec la fonction(lecture)
				Contenant tous les atomes des acides aminés du fichier pdb.
				
		Returns
		-------
			dataframe
				Contenant les cystéines et leur position pouvant faire les ponts disulfures,
					ainsi que leur distance par rapport à un cutoff (ici, 2.2 A).
	"""
	i = 0
	j = 0
	compteur = 0
	cutoff = 2.2
	liste_distance = []
	df_disulfure = pd.DataFrame((table.loc[(table['AA'] == "CYS") & (table['atom'] == "SG")]))
	disulfure = np.array(df_disulfure)
	for i in range(len(disulfure) -1):
		j = j+i
		for j in range(i+1, len(disulfure)):
			distance = (np.sqrt(((disulfure[i,5] - disulfure[(j),5])**2) + ((disulfure[i,6] \
						- disulfure[(j),6])**2) + ((disulfure[i,7] - disulfure[(j),7])**2)))
			liste_distance.append(disulfure[i,1])
			liste_distance.append(disulfure[i,0])
			liste_distance.append(disulfure[i,4])
			liste_distance.append(disulfure[j,1])
			liste_distance.append(disulfure[j,0])
			liste_distance.append(disulfure[j,4])			
			liste_distance.append(distance) 
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue', \
				 'chain', 'distance'], data = array)
	array['distance'] = array['distance'].astype(float)
	table_disulfure_cutoff = array[array["distance"] < cutoff]	
	return(table_disulfure_cutoff)

def ionic_interactions(table):
	"""
	Calcul des distances entre les atomes OE1 et OD1 des acides aminés acide aspartique et
		acide glutamique chargés négativement (COO-), et les atomes NH1, NZ et ND1 des 
		acides aminés arginine, lysine et histidine chargés positivement (NH3+).
		
		Parameters
		----------
			table : dataframe obtenu avec la fonction(lecture).
				Contenant tous les atomes des acides aminés de fichier pdb.
				
		Returns
		-------
			dataframe
				Contenant les interactions ioniques entre les acides aminés chargés opposés, 
					définit par le cutoff (ici, 6 A).
	"""
	i = 0
	j = 0
	compteur = 0
	cutoff = 6
	liste_distance = []
	df_ionic_coo = pd.DataFrame((table.loc[(((table['AA'] == "ASP") & (table['atom'] == "OD1")) | \
					(table['AA'] == "GLU") & (table['atom'] == "OE1"))]))
	ionic_coo = np.array(df_ionic_coo)
	df_ionic_nh = pd.DataFrame((table.loc[((((table['AA'] == "ARG") & (table['atom'] == "NH1")) | \
				(table['AA'] == "LYS") & (table['atom'] == "NZ")) | ((table['AA'] == "HIS") & (table['atom'] == "ND1")))]))
	ionic_nh = np.array(df_ionic_nh)
	for i in range(len(ionic_coo)):
		for j in range(len(ionic_nh)):
			distance = (np.sqrt(((ionic_coo[i,5] - ionic_nh[(j),5])**2) + ((ionic_coo[i,6] \
						- ionic_nh[(j),6])**2) + ((ionic_coo[i,7] - ionic_nh[(j),7])**2)))
			liste_distance.append(ionic_coo[i,1])
			liste_distance.append(ionic_coo[i,0])
			liste_distance.append(ionic_coo[i,4])
			liste_distance.append(ionic_nh[j,1])
			liste_distance.append(ionic_nh[j,0])
			liste_distance.append(ionic_nh[j,4])			
			liste_distance.append(distance)
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue', \
				 'chain', 'distance'], data = array)
	array['distance'] = array['distance'].astype(float)
	table_ionic_cutoff = array[array["distance"] < cutoff]
	return(table_ionic_cutoff)
		



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
    	print()
    else:
    	sys.exit("Le fichier est absent.")
    with open(Reponse,"r") as pdb_files:
    	for lines in pdb_files:
    		if lines.startswith("ATOM"):
    			comp_line=comp_line+1
    	if comp_line > 3:
    		var_lecture = lecture(Reponse)
    		disulfuric = disulphide_bridges(var_lecture)
    		print("""Intraprotein Disulphide Bridges \n Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms \n""", disulfuric)
    		print("---------------------------------------------------------------------")
    		ionic = ionic_interactions(var_lecture)
    		print("""Intraprotein Ionic Interactions \n Ionic Interactions within 6 Angstroms \n """, ionic)
    		print()
