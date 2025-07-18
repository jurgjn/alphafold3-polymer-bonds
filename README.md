# model_protein_bond

#### Pseudo code

Go through all jsons in a specified directory and parse them
Go through the bondedAtomPairs and check wether the bond is between two proteins. Collect those in a data structure "protein_protein_bonds"
+ Initialize a dictionary structure that maps every amino acid in every chain in the complex to its modified id and residue number. Initialize it with every residues own chain id and residue number. 
Go through all of the protein protein bonds, then call a model_bond_with_ligand function for each of them.

this function will check wether one of the amino acids that are forming the bond is at the end of the protein chain. So at the end or the beginnging of the chain
If it is, it will "chop off" the last amino acid and modify the current json structure with the chopped off amino acid as a ligand. 
+ The parts of the chain that are seperated will receive a new id by extending the chainId by A for the first part, by L for the ligand, and by B for the second part  
+ The dictionary structure will be updated accordingly with the new id and the modified residue number, taken into account the now shorter chain parts and the removed ligand
The bondedAtomPairs will be adjusted accordingly to model the bond now using the ligand as a "bridge". 
If it is not then it will "chop out" on of the amino acids from the bond out of its chain. then also model it as ligand but add to the bondedAtomPairs 
the bond to connect it to its original chain. And then model the bond using the ligand as a "bridge with three entrances".