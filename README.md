 # ---------------------------------------------------------------------------------------- 
				  PICI 
## Protein Interaction Calculator Implementation 
## By : Fortune Audrey & Batista Reda  
## Cursus : M2 Bioinformatique
# ----------------------------------------------------------------------------------------  


# Description 

PICI is an alternative version of a web serveur named PIC for [Protein Interaction Calculator](http://pic.mbu.iisc.ernet.in/index.html).  
PICI calculate the distance between a various number of interaction in a protein file. 
The interaction are : disulphide bridge, ionic ointeraction, main chain-main chain hydrogen bounds aromatic-aromatic interaction, aromatic-sulphur interaction, cation-n interaction. 
The program take as input a pbd file and give as output dataframe with all the interaction known in the file. 


#  How to get the script 

In order to get the program you have to :
Go on Github and click on [Here](https://github.com/AudreyyFortune/pici.git) to go on the project page. 


# How to use it  

To use the script you need : 

A python shell  
- The pdb file to be use with the hydrogen* 
- And the script of the program  
- No installation is required to use PICI. 

*After getting the pdb, you need to add the hydrogen bound to it (in PDB file hydrogen are not add in the first place).You can use Hbound, molprobity or any software of your choice.If you want to compare your results with the web serveur you must use the version without hydrogen bound. The webserveur will automaticaly in the results the hydrogen.* 


You need to have the pdb file and the script in the same folder. 
Run the program normally 
```python 
python PICI.py 
```
Then it will ask you the name of the pdb file you want to be runned.  

![](Capture4.png)

When the program as done process your pdb file the output will show you the following information : 

- The type of interaction and the cutoff. 
A dataframe with all of the principal information like : 
	- the position of the residues, the name of the residues ;
	- the chain and the distance in angstrom 

For the name of the residues, they follow the classic naming for amino acid  
To compare the results with the orginial version you can follow the link in the description adn use the web serveur.

