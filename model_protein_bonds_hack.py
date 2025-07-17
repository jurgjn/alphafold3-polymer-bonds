import json
import os
from pathlib import Path
from typing import Dict, List, Tuple, Any
from copy import deepcopy

def load_json_files(source_dir: str) -> Dict[str, Dict]:
    """
    Load all JSON files from the specified directory.
    
    Args:
        source_dir: Directory containing JSON files
        
    Returns:
        Dictionary mapping filenames to their JSON content
    """
    json_files = {}
    source_path = Path(source_dir)
    
    if not source_path.exists():
        print(f"Warning: Directory {source_dir} does not exist")
        return json_files
    
    for json_file in source_path.glob("*.json"):
        try:
            with open(json_file, 'r') as f:
                json_files[json_file.stem] = json.load(f)
            #print(f"Loaded: {json_file.name}")
        except Exception as e:
            print(f"Error loading {json_file.name}: {e}")
    
    return json_files

def initialize_residue_mapping(json_data: Dict) -> Dict[str, Dict[int, Dict[str, Any]]]:
    """
    Initialize a dictionary structure that maps every amino acid in every chain 
    to its modified id and residue number.
    
    Args:
        json_data: The AlphaFold3 JSON structure
        
    Returns:
        Dictionary mapping {chain_id: {residue_num: {modified_chain_id, modified_residue_num}}}
    """
    residue_mapping = {}
    
    for sequence in json_data.get("sequences", []):
        if "protein" in sequence:
            chain_id = sequence["protein"]["id"]
            sequence_length = len(sequence["protein"]["sequence"])
            
            residue_mapping[chain_id] = {}
            for i in range(1, sequence_length + 1):
                residue_mapping[chain_id][i] = {
                    "modified_chain_id": chain_id,
                    "modified_residue_num": i
                }
    
    return residue_mapping

def find_protein_protein_bonds(json_data: Dict) -> List[Tuple]:
    """
    Identify bonds between two protein chains from bondedAtomPairs.
    
    Args:
        json_data: The AlphaFold3 JSON structure
        
    Returns:
        List of tuples representing protein-protein bonds
    """
    protein_protein_bonds = []
    
    # Get protein chain IDs
    protein_chains = set()
    for sequence in json_data.get("sequences", []):
        if "protein" in sequence:
            protein_chains.add(sequence["protein"]["id"])
    
    # Check bondedAtomPairs for protein-protein bonds
    for bond in json_data.get("bondedAtomPairs", []):
        if len(bond) == 2:
            chain1, _, _ = bond[0]
            chain2, _, _ = bond[1]
            
            if chain1 in protein_chains and chain2 in protein_chains and chain1 != chain2:
                protein_protein_bonds.append(tuple(bond))
                #print(f"Found protein-protein bond: {chain1} -> {chain2}")
    
    return protein_protein_bonds

def get_sequence_info(json_data: Dict, chain_id: str) -> Dict:
    """
    Get sequence information for a specific chain.
    
    Args:
        json_data: The AlphaFold3 JSON structure
        chain_id: Chain identifier
        
    Returns:
        Dictionary with sequence information
    """
    for sequence in json_data.get("sequences", []):
        if "protein" in sequence and sequence["protein"]["id"] == chain_id:
            return sequence["protein"]
    return {}

def is_terminal_residue(sequence: str, residue_position: int) -> bool:
    """
    Check if a residue is at the terminal ends of a protein chain.
    
    Args:
        sequence: Protein sequence
        residue_position: 1-based position of the residue
        
    Returns:
        True if residue is at N-terminus (position 1) or C-terminus (last position)
    """
    return residue_position == 1 or residue_position == len(sequence)

def update_residue_mapping_for_terminal_split(residue_mapping: Dict, chain_id: str, 
                                            split_position: int, sequence_length: int,
                                            is_c_terminal: bool) -> None:
    """
    Update residue mapping when splitting a chain at terminal position.
    
    Args:
        residue_mapping: The residue mapping dictionary to update
        chain_id: Original chain ID
        split_position: Position where the split occurs
        sequence_length: Length of the original sequence
        is_c_terminal: True if splitting at C-terminal, False for N-terminal
    """
    if is_c_terminal:
        # C-terminal split: chain keeps positions 1 to split_position-1
        # Ligand gets the last residue
        for pos in range(1, split_position):
            residue_mapping[chain_id][pos]["modified_chain_id"] = f"{chain_id}A"
            residue_mapping[chain_id][pos]["modified_residue_num"] = pos
        
        # The terminal residue becomes ligand
        residue_mapping[chain_id][split_position]["modified_chain_id"] = f"{chain_id}L"
        residue_mapping[chain_id][split_position]["modified_residue_num"] = 1
    else:
        # N-terminal split: first residue becomes ligand, rest shifts down
        residue_mapping[chain_id][1]["modified_chain_id"] = f"{chain_id}L"
        residue_mapping[chain_id][1]["modified_residue_num"] = 1
        
        for pos in range(2, sequence_length + 1):
            residue_mapping[chain_id][pos]["modified_chain_id"] = f"{chain_id}A"
            residue_mapping[chain_id][pos]["modified_residue_num"] = pos - 1

def update_residue_mapping_for_internal_split(residue_mapping: Dict, chain_id: str, 
                                            split_position: int, sequence_length: int) -> None:
    """
    Update residue mapping when splitting a chain at internal position.
    
    Args:
        residue_mapping: The residue mapping dictionary to update
        chain_id: Original chain ID
        split_position: Position where the split occurs
        sequence_length: Length of the original sequence
    """
    # Part A: residues 1 to split_position-1
    for pos in range(1, split_position):
        residue_mapping[chain_id][pos]["modified_chain_id"] = f"{chain_id}A"
        residue_mapping[chain_id][pos]["modified_residue_num"] = pos
    
    # Ligand: the residue at split_position
    residue_mapping[chain_id][split_position]["modified_chain_id"] = f"{chain_id}L"
    residue_mapping[chain_id][split_position]["modified_residue_num"] = 1
    
    # Part B: residues split_position+1 to end
    for pos in range(split_position + 1, sequence_length + 1):
        residue_mapping[chain_id][pos]["modified_chain_id"] = f"{chain_id}B"
        residue_mapping[chain_id][pos]["modified_residue_num"] = pos - split_position

def model_bond_with_ligand(json_data: Dict, bond: Tuple, residue_mapping: Dict) -> Dict:
    """
    Modify JSON structure to model protein-protein bond using a ligand bridge.
    
    Args:
        json_data: Original AlphaFold3 JSON structure
        bond: Tuple representing the protein-protein bond
        residue_mapping: Dictionary tracking residue mappings
        
    Returns:
        Modified JSON structure with ligand bridge
    """
    # Deep copy to avoid modifying original
    modified_json = deepcopy(json_data)
    
    # Extract bond information
    atom1, atom2 = bond
    chain1_id, seq_num1, atom_name1 = atom1
    chain2_id, seq_num2, atom_name2 = atom2
    
    # Get sequence information
    chain1_info = get_sequence_info(json_data, chain1_id)
    chain2_info = get_sequence_info(json_data, chain2_id)
    
    if not chain1_info or not chain2_info:
        print(f"Warning: Could not find sequence info for chains {chain1_id} or {chain2_id}")
        return modified_json
    
    # Determine which residue to use as ligand (prefer terminal residues)
    chain1_sequence = chain1_info["sequence"]
    chain2_sequence = chain2_info["sequence"]
    
    chain1_is_terminal = is_terminal_residue(chain1_sequence, seq_num1)
    chain2_is_terminal = is_terminal_residue(chain2_sequence, seq_num2)
    
    # Prefer terminal residues, if both or neither are terminal, use chain1
    use_chain1_as_ligand = chain1_is_terminal or not chain2_is_terminal
    
    # Convert single letter amino acid to CCD code
    aa_to_ccd = {
        'G': 'GLY', 'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE',
        'P': 'PRO', 'F': 'PHE', 'Y': 'TYR', 'W': 'TRP', 'S': 'SER',
        'T': 'THR', 'C': 'CYS', 'M': 'MET', 'N': 'ASN', 'Q': 'GLN',
        'D': 'ASP', 'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'
    }
    
    new_sequences = []
    new_bonded_pairs = []
    
    if use_chain1_as_ligand:
        ligand_chain_id = chain1_id
        ligand_seq_num = seq_num1
        ligand_residue = chain1_sequence[seq_num1 - 1]
        ligand_ccd = aa_to_ccd.get(ligand_residue, 'UNK')
        
        # Update residue mapping
        if chain1_is_terminal:
            is_c_terminal = seq_num1 == len(chain1_sequence)
            update_residue_mapping_for_terminal_split(residue_mapping, chain1_id, 
                                                    seq_num1, len(chain1_sequence), is_c_terminal)
            
            # Create modified sequence
            if is_c_terminal:
                modified_sequence = chain1_sequence[:-1]
                new_chain_id = f"{chain1_id}A"
            else:
                modified_sequence = chain1_sequence[1:]
                new_chain_id = f"{chain1_id}A"
        else:
            # Internal residue
            update_residue_mapping_for_internal_split(residue_mapping, chain1_id, 
                                                    seq_num1, len(chain1_sequence))
            
            # Split into two parts
            part_a_sequence = chain1_sequence[:seq_num1-1]
            part_b_sequence = chain1_sequence[seq_num1:]
            
            # Add part A if it exists
            if part_a_sequence:
                new_sequences.append({
                    "protein": {
                        "id": f"{chain1_id}A",
                        "sequence": part_a_sequence
                    }
                })
                # Bond from part A to ligand
                new_bonded_pairs.append([[f"{chain1_id}A", len(part_a_sequence), "C"], 
                                       [f"{chain1_id}L", 1, "N"]])
            
            # Add part B if it exists
            if part_b_sequence:
                new_sequences.append({
                    "protein": {
                        "id": f"{chain1_id}B",
                        "sequence": part_b_sequence
                    }
                })
                # Bond from ligand to part B
                new_bonded_pairs.append([[f"{chain1_id}L", 1, "C"], 
                                       [f"{chain1_id}B", 1, "N"]])
        
        # Add ligand
        ligand_id = f"{chain1_id}L"
        new_sequences.append({
            "ligand": {
                "id": ligand_id,
                "ccdCodes": [ligand_ccd]
            }
        })
        
        # For terminal case, add the remaining protein sequence
        if chain1_is_terminal and modified_sequence:
            new_sequences.append({
                "protein": {
                    "id": new_chain_id,
                    "sequence": modified_sequence
                }
            })
            
            # Add bonds for terminal case
            if seq_num1 == len(chain1_sequence):  # C-terminal
                new_bonded_pairs.append([[new_chain_id, len(modified_sequence), "C"], 
                                       [ligand_id, 1, "N"]])
            else:  # N-terminal
                new_bonded_pairs.append([[ligand_id, 1, "C"], 
                                       [new_chain_id, 1, "N"]])
        
        # Add other sequences unchanged
        for sequence in modified_json["sequences"]:
            if "protein" in sequence and sequence["protein"]["id"] != chain1_id:
                new_sequences.append(sequence)
            elif "ligand" in sequence or "dna" in sequence or "rna" in sequence:
                new_sequences.append(sequence)
        
        # Bond from ligand to chain2 (using original mapping for chain2)
        chain2_mapped = residue_mapping[chain2_id][seq_num2]
        new_bonded_pairs.append([[ligand_id, 1, atom_name1], 
                               [chain2_mapped["modified_chain_id"], 
                                chain2_mapped["modified_residue_num"], atom_name2]])
    
    # Add existing bonds (excluding the original protein-protein bond)
    for existing_bond in modified_json.get("bondedAtomPairs", []):
        if tuple(existing_bond) != bond:
            new_bonded_pairs.append(existing_bond)
    
    modified_json["sequences"] = new_sequences
    modified_json["bondedAtomPairs"] = new_bonded_pairs
    
    return modified_json

def process_json_files(source_dir: str, output_dir: str = "output/jsons/ubn_links_modified/") -> None:
    """
    Process all JSON files in the source directory and create modified versions.
    
    Args:
        source_dir: Directory containing original JSON files
        output_dir: Directory to save modified JSON files
    """
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    # Load all JSON files
    json_files = load_json_files(source_dir)
    
    for filename, json_data in json_files.items():
        #print(f"\nProcessing {filename}...")
        
        # Initialize residue mapping
        residue_mapping = initialize_residue_mapping(json_data)
        
        # Find protein-protein bonds
        protein_bonds = find_protein_protein_bonds(json_data)
        
        if not protein_bonds:
            print(f"No protein-protein bonds found in {filename}")
            continue
        
        # Process each protein-protein bond
        modified_json = json_data
        for bond in protein_bonds:
            #print(f"Modifying bond: {bond}")
            modified_json = model_bond_with_ligand(modified_json, bond, residue_mapping)
        
        # Save modified JSON
        output_path = Path(output_dir) / f"{filename}_modified.json"
        with open(output_path, 'w') as f:
            json.dump(modified_json, f, indent=2)
        
        #print(f"Saved modified file: {output_path}")

# Set source directory
source_dir = "output/jsons/ubn_links/"

# Process all JSON files
process_json_files(source_dir)