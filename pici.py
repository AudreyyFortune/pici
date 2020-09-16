""" 
Sujet 6 : Protein Interactions Calculator implémentation.

----------------------------------------------------------------
BATISTA Reda &  FORTUNE Audrey
Cursus : M2 Bio informatique
----------------------------------------------------------------
"""

def lecture(name):
	"""
	Reading a pdb file and extracting the data.

		Parameters 
		---------- 
		name : pdb file
    		Contain the atomic coordinates of the amino acid atoms.

		Returns
		------- 
		dataframe 
    		Contain all the lines of the pdb file containing the amino acids.
    		Table with as columns: the residue, the position, the atom, the position of 
    		the atom, the chain, the coordinates x, y and z, as well as the occupancy.
	"""
	flag_count = 0
	flag_shape = 0
	with open(name, "r") as pdb_file:
		for line in pdb_file:
			if line.startswith("ATOM"):
				flag_shape += 1
	p = np.full((flag_shape,9),np.nan)
	table = pd.DataFrame(p)
	table.columns = ['AA','position_AA','atom','position_atom','chain','coord_x','coord_y',
					'coord_z','occupancy']
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
				table.iloc[flag_count] = [AA_name,AA_position,atom_name,atom_position,
											chain,x,y,z,occupancy]
				flag_count += 1
				table['position_atom'] = table['position_atom'].replace(['0'],np.nan)
		return(table) 


def disulphide_bridges(table):
	"""
	Calculation of the distances between each SG atom of the cysteines, and selection of 
		the disulfide bridges.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
				
		Returns
		-------
			dataframe
				Contain cysteines and their position that can make the disulfide bridges, 
				as well as their distance from a cutoff (here, 2.2 A).
	"""
	i = 0
	j = 0
	compteur = 0
	cutoff = 2.2
	liste_distance = []
	df_disulfure = pd.DataFrame((table.loc[(table['AA'] == "CYS") & (table['atom'] == 
					"SG")]))
	disulfure = np.array(df_disulfure)
	for i in range(len(disulfure) -1):
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
	array = array.reshape((compteur, 7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue',
				 'chain', 'distance'], data = array)
	array['distance'] = array['distance'].astype(float)
	array['distance'] = round(array['distance'], 2)
	table_disulfure_cutoff = array[array["distance"] < cutoff]	
	return(table_disulfure_cutoff)


def ionic_interaction(table):
	"""
	Calculation of the distances between the OE1 or OD1 atoms of the amino acids aspartic 
		acid or glutamic acid negatively charge (COO-), and the NH1, NZ or ND1 atoms of the 
		positively charged amino acids arginine, lysine or histidine (NH3+).
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
				
		Returns
		-------
			dataframe
				Contain the ionic interactions between oppositely charged amino acids, 
				defined by the cutoff (here, 6 A).
	"""
	i = 0
	j = 0
	compteur = 0
	cutoff = 6
	liste_distance = []
	df_ionic_coo = pd.DataFrame((table.loc[(((table['AA'] == "ASP") & (table['atom'] == 
					"OD1")) | (table['AA'] == "GLU") & (table['atom'] == "OE1"))]))
	ionic_coo = np.array(df_ionic_coo)
	df_ionic_nh = pd.DataFrame((table.loc[((((table['AA'] == "ARG") & (table['atom'] == 
					"NH1")) | (table['AA'] == "LYS") & (table['atom'] == "NZ")) | 
					((table['AA'] == "HIS") & (table['atom'] == "ND1")))]))
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
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue',
				 'chain', 'distance'], data = array)
	array['distance'] = array['distance'].astype(float)
	array['distance'] = round(array['distance'], 2)
	table_ionic_cutoff = array[array["distance"] < cutoff]
	return(table_ionic_cutoff)
		

def mc_mc_interaction(table):
	"""
	Calculate the distances of the hydrogen interactions between the mains tailings chains.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
				
		Returns
		-------
			dataframe
				Contain the interactions between the main residue chains, defined by the 
				cutoff (here, 3.5 A).
	"""
	i, j, k, l = 0, 0, 0, 0
	compteur = 0
	cutoff = 3.5
	liste_distance_no = []
	liste_distance_ho = []
	df_donor = pd.DataFrame((table.loc[(table['atom'] == "N")]))
	df_donor = df_donor.iloc[::-1]
	donor = np.array(df_donor)
	df_acceptor = pd.DataFrame((table.loc[(table['atom'] == "O")]))
	acceptor = np.array(df_acceptor)
	df_donor_h = pd.DataFrame((table.loc[(table['atom'] == "H")]))
	df_donor_h = df_donor_h.iloc[::-1]
	donor_h = np.array(df_donor_h)
	for i in range(len(donor))	:
		for j in range(len(acceptor) -(i+2)):
			distance_no = (np.sqrt(((donor[i,5] - acceptor[(j),5])**2) + ((donor[i,6] \
						- acceptor[(j),6])**2) + ((donor[i,7] - acceptor[(j),7])**2)))
			liste_distance_no.append(donor[i,1])
			liste_distance_no.append(donor[i,4])
			liste_distance_no.append(donor[i,0])
			liste_distance_no.append(donor[i,2])
			liste_distance_no.append(acceptor[j,1])
			liste_distance_no.append(acceptor[j,4])
			liste_distance_no.append(acceptor[j,0])	
			liste_distance_no.append(acceptor[j,2])
			liste_distance_no.append(distance_no)
			compteur += 1
	for k in range(len(donor_h)):
		for l in range(len(acceptor) -(k+2)):
			distance_ho = (np.sqrt(((donor_h[k,5] - acceptor[(l),5])**2) + ((donor_h[k,6] \
						- acceptor[(l),6])**2) + ((donor_h[k,7] - acceptor[(l),7])**2)))
			liste_distance_ho.append(distance_ho)
	array = np.array(liste_distance_no)
	array = array.reshape((compteur,9))
	df = pd.DataFrame(columns = ['position', 'chain', 'residu', 'atom', 'position', 
							'chain', 'residu', 'atom', 'Dd-a'], data = array)
	df['Dd-a'] = df['Dd-a'].astype(float)
	df['Dd-a'] = round(df['Dd-a'], 2)
	df_ho = pd.DataFrame(liste_distance_ho, columns = ['Dh-a'])
	df_ho['Dh-a'] = round(df_ho['Dh-a'], 2)
	df_no_ho = pd.concat([df,df_ho], axis = 1)
	table_mc_mc_cutoff = df_no_ho[df_no_ho["Dd-a"] < cutoff]
	return(table_mc_mc_cutoff)
	

def moy_coord_1c(table):
	"""
	Averaging the x, y and z coordinates of the six ring atoms of the aromatic amino acids
		tyrosine and phenylalanine.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
					Containing all the amino acid atoms in the pdb file.
				
		Returns
		-------
			list
				Contain the means of the coordinates of the aromatic rings of tyrosine 
				and phenylalanine.
	"""
	max_1c = 6
	dataframe_1c = []
	somme, element = 0, 0
	liste_moy_1c = []
	df_1cycle = pd.DataFrame((table.loc[((table['AA'] == "TYR") | (table['AA'] == "PHE"))
				 & ((table['atom'] == "CD1") | (table['atom'] == "CD2") | (table['atom'] 
				== "CE1") | (table['atom'] == "CZ") | (table['atom'] == "CE2") | 
				(table['atom'] == "CG") | (table['atom'] == "NE1") | (table['atom'] == 
				"CZ2") | (table['atom'] == "CH2") | (table['atom'] == "CZ3") | 
				(table['atom'] == "CE3"))]))
	df_1c = df_1cycle[['coord_x','coord_y','coord_z']]
	array_1c = np.array(df_1c)
	x_1c = array_1c[:,0]
	y_1c = array_1c[:,1]
	z_1c = array_1c[:,2]
	while len(x_1c) > max_1c:
		top_x = x_1c[:max_1c]
		dataframe_1c.append(top_x)
		x_1c = x_1c[max_1c:]
	else:
		dataframe_1c.append(x_1c)
		while len(y_1c) > max_1c:
			top_y = y_1c[:max_1c]
			dataframe_1c.append(top_y)
			y_1c = y_1c[max_1c:]
		else:
			dataframe_1c.append(y_1c)
			while len(z_1c) > max_1c:
				top_z = z_1c[:max_1c]
				dataframe_1c.append(top_z)
				z_1c = z_1c[max_1c:]
			else:
				dataframe_1c.append(z_1c)
	if len(dataframe_1c) > 3 : 
		for elements in range(len(dataframe_1c)):
			moy = dataframe_1c[element].mean()
			liste_moy_1c.append(moy)
			element += 1
	return(liste_moy_1c)
				
def moy_coord_2c(table):
	"""
	Averaging the x, y and z coordinates of the nine ring atoms of the aromatic amino 
		acids tryptophan.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
					Containing all the amino acid atoms in the pdb file.
				
		Returns
		-------
			list
				Contain the mean of the coordinates of the aromatic rings of tryptophan.
	"""
	max_2c = 9
	dataframe_2c = []
	somme, element = 0, 0
	liste_moy_2c = []
	df_2cycles = pd.DataFrame((table.loc[(table['AA'] == "TRP") & ((table['atom'] == 
				"CD1") | (table['atom'] == "CD2") | (table['atom'] == "CE1") | 
				(table['atom'] == "CZ") | (table['atom'] == "CE2") | (table['atom'] == 
				"CG") | (table['atom'] == "NE1") | (table['atom'] == "CZ2") | 
				(table['atom'] == "CH2") | (table['atom'] == "CZ3") | (table['atom'] == 
				"CE3"))]))
	df_2c = df_2cycles[['coord_x','coord_y','coord_z']]
	array_2c = np.array(df_2c)
	x_2c = array_2c[:,0]
	y_2c = array_2c[:,1]
	z_2c = array_2c[:,2]
	while len(x_2c) > max_2c:
		top_x = x_2c[:max_2c]
		dataframe_2c.append(top_x)
		x_2c = x_2c[max_2c:]
	else:
		dataframe_2c.append(x_2c)
		while len(y_2c) > max_2c:
			top_y = y_2c[:max_2c]
			dataframe_2c.append(top_y)
			y_2c = y_2c[max_2c:]
		else:
			dataframe_2c.append(y_2c)
			while len(z_2c) > max_2c:
				top_z = z_2c[:max_2c]
				dataframe_2c.append(top_z)
				z_2c = z_2c[max_2c:]
			else:
				dataframe_2c.append(z_2c)
	if len(dataframe_2c) > 3 : 
		for elements in range(len(dataframe_2c)):
			moy = dataframe_2c[element].mean()
			liste_moy_2c.append(moy)
			element += 1
	return(liste_moy_2c)


def aromatic_interaction(table, liste_moy_1c,liste_moy_2c):
	"""
	Calculation of the distances between the average of the atoms of the aromatic rings 
		and selection of the aromatic interactions.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
			liste_moy_1c: list
				Containing the means of the coordinates of the aromatic rings of tyrosine 
				and phenylalanine.
			liste_moy_2c: list
				Contain the mean of the coordinates of the aromatic rings of tryptophan.
				
		Returns
		-------
			dataframe
				Contain the interactions between aromatic amino acids, defined by a cutoff
				between 4.5 and 7 A.
	"""
	i, j = 0, 0
	compteur = 0
	cutoff_min = 4.5
	cutoff_max = 7
	liste_distance = []
	df_aromatic = pd.DataFrame((table.loc[((table['AA'] == "TYR") | (table['AA'] == "PHE")
				 | (table['AA'] == "TRP")) & (table['atom'] == "CG")]))
	df_aromatic = df_aromatic[['AA','position_AA','chain']]
	aromatic = np.array(df_aromatic)
	liste_1c = np.array(liste_moy_1c)
	liste_1c.resize(3,len(aromatic))
	liste_1c = np.transpose(liste_1c)
	aromatic = np.concatenate((aromatic,liste_1c), axis=1)
	for i in range(len(aromatic) -1):
		for j in range(i+1, len(aromatic)):
			distance = (np.sqrt(((aromatic[i,3] - aromatic[(j),3])**2) + ((aromatic[i,4] \
						- aromatic[(j),4])**2) + ((aromatic[i,5] - aromatic[(j),5])**2)))
			liste_distance.append(aromatic[i,1])
			liste_distance.append(aromatic[i,0])
			liste_distance.append(aromatic[i,2])
			liste_distance.append(aromatic[j,1])
			liste_distance.append(aromatic[j,0])
			liste_distance.append(aromatic[j,2])			
			liste_distance.append(distance) 
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue',
				 'chain', 'D(centroid-centroid)'], data = array)
	array['D(centroid-centroid)'] = array['D(centroid-centroid)'].astype(float)
	array['D(centroid-centroid)'] = round(array['D(centroid-centroid)'], 2)
	table_aromatic_cutoff = array[ (array["D(centroid-centroid)"] > cutoff_min) & \
							(array["D(centroid-centroid)"] < cutoff_max)]	
	return(table_aromatic_cutoff)


def aromatic_sulfur_interaction(table, liste_moy_1c, liste_moy_2c):
	"""
	Calculation of distances between the average of the atoms of the aromatic rings and 
		the atoms sulphur of cysteine and methionine.
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
			liste_moy_1c: liste
				Containing the means of the coordinates of the aromatic rings of tyrosine 
				and phenylalanine.
			liste_moy_2c: liste
				Contain the mean of the coordinates of the aromatic rings of tryptophan.
				
		Returns
		-------
			dataframe
				Contain the interactions between aromatic acids amino and sulphide amino 
				acids, defined by a cutoff of 5.3 A.
	"""
	i, j = 0, 0
	compteur = 0
	cutoff = 5.3
	liste_distance = []
	df_sulfur = pd.DataFrame((table.loc[((table['AA'] == "CYS") & (table['atom'] == "SG"))
					| ((table['AA'] == "MET") & (table['atom'] == "SD"))]))
	sulfur = np.array(df_sulfur)
	df_aromatic = pd.DataFrame((table.loc[((table['AA'] == "TYR") | (table['AA'] == "PHE")
				 | (table['AA'] == "TRP")) & (table['atom'] == "CG")]))
	df_aromatic = df_aromatic[['AA','position_AA','chain']]
	aromatic = np.array(df_aromatic)
	aromatic = np.array(df_aromatic)
	liste_1c = np.array(liste_moy_1c)
	liste_1c.resize(3,len(aromatic))
	liste_1c = np.transpose(liste_1c)
	aromatic = np.concatenate((aromatic,liste_1c), axis=1)
	for i in range(len(sulfur)):
		for j in range(len(aromatic)):
			distance = (np.sqrt(((sulfur[i,5] - aromatic[(j),3])**2) + ((sulfur[i,6] \
						- aromatic[(j),4])**2) + ((sulfur[i,7] - aromatic[(j),5])**2)))
			liste_distance.append(aromatic[j,1])
			liste_distance.append(aromatic[j,0])
			liste_distance.append(aromatic[j,2])
			liste_distance.append(sulfur[i,1])
			liste_distance.append(sulfur[i,0])
			liste_distance.append(sulfur[i,4])			
			liste_distance.append(distance)
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue',
				 'chain', 'D(centroid-sulfur)'], data = array)
	array['D(centroid-sulfur)'] = array['D(centroid-sulfur)'].astype(float)
	array['D(centroid-sulfur)'] = round(array['D(centroid-sulfur)'], 2)
	table_aromatic_sulfur_cutoff = array[array["D(centroid-sulfur)"] < cutoff]
	return(table_aromatic_sulfur_cutoff)
	
	
def cation_pi_interaction(table, liste_moy_1c, liste_moy_2c):
	"""
	Calculation of the distances between the average of the atoms of the aromatic rings 
		and the cationic residues (arginine and lysine).
		
		Parameters
		----------
			table : dataframe obtained with the function(lecture)
				Containing all the amino acid atoms in the pdb file.
			liste_moy_1c: liste
				Containing the means of the coordinates of the aromatic rings of tyrosine 
				and phenylalanine.
			liste_moy_2c: liste
				Contain the mean of the coordinates of the aromatic rings of tryptophan.
				
		Returns
		-------
			dataframe
				Contain the interactions between aromatic amino acids and cationic amino
				acids, defined by a 6 A cutoff.
	"""
	i = 0
	j = 0
	compteur = 0
	cutoff = 6
	liste_distance = []
	df_ionic = pd.DataFrame((table.loc[((table['AA'] == "ARG") & (table['atom'] == "NH1")) 
				| ((table['AA'] == "LYS") & (table['atom'] == "NZ"))]))
	ionic = np.array(df_ionic)
	df_aromatic = pd.DataFrame((table.loc[((table['AA'] == "TYR") | (table['AA'] == "PHE")
				 | (table['AA'] == "TRP")) & (table['atom'] == "CG")]))
	df_aromatic = df_aromatic[['AA','position_AA','chain']]
	aromatic = np.array(df_aromatic)
	aromatic = np.array(df_aromatic)
	liste_1c = np.array(liste_moy_1c)
	liste_1c.resize(3,len(aromatic))
	liste_1c = np.transpose(liste_1c)
	aromatic = np.concatenate((aromatic,liste_1c), axis=1)
	for i in range(len(ionic)):
		for j in range(len(aromatic)):
			distance = (np.sqrt(((ionic[i,5] - aromatic[(j),3])**2) + ((ionic[i,6] \
						- aromatic[(j),4])**2) + ((ionic[i,7] - aromatic[(j),5])**2)))
			liste_distance.append(aromatic[j,1])
			liste_distance.append(aromatic[j,0])
			liste_distance.append(aromatic[j,2])
			liste_distance.append(ionic[i,1])
			liste_distance.append(ionic[i,0])
			liste_distance.append(ionic[i,4])			
			liste_distance.append(distance)
			compteur += 1
	array = np.array(liste_distance)
	array = array.reshape((compteur,7))
	array = pd.DataFrame(columns = ['position', 'residue', 'chain', 'position', 'residue',
				 'chain', 'D(cation-pi)'], data = array)
	array['D(cation-pi)'] = array['D(cation-pi)'].astype(float)
	array['D(cation-pi)'] = round(array['D(cation-pi)'], 2)
	table_cation_pi_cutoff = array[array["D(cation-pi)"] < cutoff]
	return(table_cation_pi_cutoff)


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
    		print("---------------------------------------------------------------------")
    		print("Intraprotein Disulphide Bridges \n ")
    		print("Disulphide bridges: Between sulphur atoms of cysteines within 2.2 Angstroms \n",
    				 disulfuric)
    		print("---------------------------------------------------------------------")
    		print()
    		ionic = ionic_interaction(var_lecture)
    		print("Intraprotein Ionic Interactions \n")
    		print("Ionic Interactions within 6 Angstroms \n ", ionic)
    		print("---------------------------------------------------------------------")
    		print()
    		mc_mc = mc_mc_interaction(var_lecture)
    		print("Intraprotein Main Chain-Main Chain Hydrogen Bonds \n")
    		print(" 	DONOR 			ACCEPTOR 		PARAMETERS \n" ,mc_mc)
    		print("---------------------------------------------------------------------")
    		print()
    		moy_coordonne_1c = moy_coord_1c(var_lecture)
    		moy_coordonne_2c = moy_coord_2c(var_lecture)
    		aromatic = aromatic_interaction(var_lecture,moy_coordonne_1c,moy_coordonne_2c)
    		print("Intraprotein Aromatic-Aromatic Interactions \n")
    		print("Aromatic-Aromatic Interactions within 4.5 and 7 Angstroms \n",aromatic)
    		print("---------------------------------------------------------------------")
    		print()
    		aromatic_sulfur = aromatic_sulfur_interaction(var_lecture,moy_coordonne_1c,moy_coordonne_2c)
    		print("Intraprotein Aromatic-Sulphur Interactions  \n")
    		print("Aromatic-Sulphur Interactions within 5.3 Angstroms \n",aromatic_sulfur)
    		print("---------------------------------------------------------------------")
    		print()
    		cation_pi = cation_pi_interaction(var_lecture,moy_coordonne_1c,moy_coordonne_2c)
    		print("Intraprotein Cation-Pi Interactions \n")
    		print("Cation-Pi Interactions within 6 Angstroms \n",cation_pi)
    		print()
    		
    		
    		
    		
