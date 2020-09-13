def lecture(name):
    flag_count=0
    flag_shape=0 
    with open(name, "r") as pdb_file: 
        for line in pdb_file:
            if line.startswith("ATOM"): 
                flag_shape+=1   
    p = np.full((flag_shape,9),np.nan) 
    table = pd.DataFrame(p)
    table.columns=["AA","position_AA","atome","position_atome","chaine","coordx","coordy","coordz","occupancy"]    
    with open(name, "r") as pdb_file: 
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
                table.iloc[flag_count]=[AA_name,AA_position,atome_name,atome_position,chaine,x,y,z,occupancy] 
                flag_count+=1
        return(table) 

def creation_table(table): 
    table_souffre=pd.DataFrame((table.loc[(table['AA']=="CYS") | (table['AA'] =="MET")])) 
    table_hydrophobic=pd.DataFrame((table.loc[(table['AA']=="ALA") | (table['AA'] =="VAL") | (table['AA'] =="LEU") | (table['AA'] =="ILE") | (table['AA'] =="MET") | (table['AA'] =="PHE") | (table['AA'] =="TRP") | (table['AA'] =="PRO") | (table['AA'] =="TYR")])) 
    
    



if __name__ == "__main__": 
    import numpy as np
    import pandas as pd 
    var_lecture= lecture("2gkt.pdb") 
    var_creation_table = creation_table(var_lecture)
