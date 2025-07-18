# model_protein_bond

## Pseudo Code

### 1. Parse Input Files (`load_json_files()`)
- Go through all JSONs in a specified directory and parse them
- Go through the `bondedAtomPairs` and check whether the bond is between two proteins (`find_protein_protein_bonds()`)
- Collect protein-protein bonds in a data structure called `protein_protein_bonds`

### 2. Initialize Data Structure (`initialize_residue_mapping()`)
- Initialize a dictionary that maps every amino acid in every chain to its modified ID and residue number
- Initialize with each residue's own chain ID and residue number

### 3. Process Protein-Protein Bonds (`model_bond_with_ligand()`)
For each protein-protein bond, call `model_bond_with_ligand()` function which calls `process_chain_bond()`:

#### Case A: Terminal Amino Acid Bond (`is_terminal_residue()`)
If one amino acid is at the end/beginning of the protein chain:
- **Chop off** the terminal amino acid and model it as a ligand (`create_ligand_from_residue()`)
- **Update chain IDs** (`create_protein_sequence()`): 
  - First part: `chainId + A`
  - Ligand: `chainId + L` 
  - Second part: `chainId + B`
- **Update dictionary** (`update_residue_mapping_for_terminal_split()`) with new IDs and adjusted residue numbers
- **Adjust bondedAtomPairs** (`correct_chain_and_resnum()`, `add_peptide_bond()`) to model the bond using the ligand as a bridge

#### Case B: Internal Amino Acid Bond  
If the amino acid is internal to the chain:
- **Chop out** the amino acid and model it as a ligand (`create_ligand_from_residue()`)
- **Add bonds** (`add_peptide_bond()`) to connect the ligand to its original chain
- **Update mapping** (`update_residue_mapping_for_internal_split()`)
- **Model the bond** using the ligand as a "bridge with three connections"

### 4. Main Processing (`process_json_files()`, `main()`)
- Orchestrates the entire workflow for all JSON files in the input directory