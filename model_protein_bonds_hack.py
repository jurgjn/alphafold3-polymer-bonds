#!/usr/bin/env python3
"""
Protein Bond Modeling Script

This script processes AlphaFold3 JSON files to model protein-protein bonds 
using ligand bridges. It identifies bonds between protein chains and modifies 
the structure to represent these bonds through intermediate ligand molecules.

Usage:
    python model_protein_bonds_hack.py --source-dir input/ --output-dir output/
    python model_protein_bonds_hack.py -s input/ -o output/ --verbose

Author: Generated from Jupyter notebook
"""

import json
import os
import argparse
import traceback
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
            print(f"Loaded: {json_file.name}")
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
    protein_bonds = []
    
    # Get protein chain IDs
    protein_chains = set()
    for sequence in json_data.get("sequences", []):
        if "protein" in sequence:
            protein_chains.add(sequence["protein"]["id"])
    
    # Check bondedAtomPairs for protein bonds (both inter and intra-chain)
    for bond in json_data.get("bondedAtomPairs", []):
        if len(bond) == 2:
            chain1, _, _ = bond[0]
            chain2, _, _ = bond[1]
            
            # Include bonds between protein chains AND within the same protein chain
            if chain1 in protein_chains and chain2 in protein_chains:
                protein_bonds.append(tuple(bond))
                if chain1 == chain2:
                    print(f"Found intra-chain bond: {chain1} (internal)")
                else:
                    print(f"Found inter-chain bond: {chain1} -> {chain2}")
    
    return protein_bonds

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

def update_residue_mapping_for_terminal_split(residue_mapping: Dict, chain_id_orig: str, 
                                            split_position: int, sequence_length: int,
                                            is_c_terminal: bool) -> None:
    """
    Update residue mapping when splitting a chain at terminal position.
    
    Args:
        residue_mapping: The residue mapping dictionary to update
        chain_id_orig: Original chain ID
        split_position: Position where the split occurs
        sequence_length: Length of the original sequence
        is_c_terminal: True if splitting at C-terminal, False for N-terminal
    """
    if is_c_terminal:
        # C-terminal split: chain keeps positions 1 to split_position-1
        # Ligand gets the last residue
        for pos in range(1, split_position):
            residue_mapping[chain_id_orig][pos]["modified_chain_id"] = f"{chain_id_orig}A"
            residue_mapping[chain_id_orig][pos]["modified_residue_num"] = pos
        
        # The terminal residue becomes ligand
        residue_mapping[chain_id_orig][split_position]["modified_chain_id"] = f"{chain_id_orig}L"
        residue_mapping[chain_id_orig][split_position]["modified_residue_num"] = 1
    else:
        # N-terminal split: first residue becomes ligand, rest shifts down
        residue_mapping[chain_id_orig][1]["modified_chain_id"] = f"{chain_id_orig}L"
        residue_mapping[chain_id_orig][1]["modified_residue_num"] = 1
        
        for pos in range(2, sequence_length + 1):
            residue_mapping[chain_id_orig][pos]["modified_chain_id"] = f"{chain_id_orig}A"
            residue_mapping[chain_id_orig][pos]["modified_residue_num"] = pos - 1

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
    Modify JSON structure to model protein bond using a ligand bridge.
    Handles both inter-chain and intra-chain bonds.
    
    Args:
        json_data: Original AlphaFold3 JSON structure
        bond: Tuple representing the protein bond
        residue_mapping: Dictionary tracking residue mappings
        
    Returns:
        Modified JSON structure with ligand bridge
    """
    # Deep copy to avoid modifying original
    modified_json = deepcopy(json_data)
    
    # Extract bond information
    atom1, atom2 = bond
    chain1_id_orig, seq_num1_old, atom_name1 = atom1
    chain2_id_orig, seq_num2_old, atom_name2 = atom2
    
    # Check if this is an intra-chain bond
    is_intra_chain = chain1_id_orig == chain2_id_orig
    
    if is_intra_chain:
        print(f"Processing intra-chain bond in {chain1_id_orig}: residue {seq_num1_old} to {seq_num2_old}")
    else:
        print(f"Processing inter-chain bond: {chain1_id_orig} to {chain2_id_orig}")
    return process_chain_bond(modified_json, bond, is_intra_chain, residue_mapping, json_data)

def get_amino_acid_ccd_map():
    """Return the amino acid to CCD code mapping."""
    return {
        'G': 'GLY', 'A': 'ALA', 'V': 'VAL', 'L': 'LEU', 'I': 'ILE',
        'P': 'PRO', 'F': 'PHE', 'Y': 'TYR', 'W': 'TRP', 'S': 'SER',
        'T': 'THR', 'C': 'CYS', 'M': 'MET', 'N': 'ASN', 'Q': 'GLN',
        'D': 'ASP', 'E': 'GLU', 'K': 'LYS', 'R': 'ARG', 'H': 'HIS'
    }

def create_ligand_from_residue(chain_id: str, residue_char: str, ligand_suffix: str = "L") -> Dict:
    """Create a ligand sequence entry from an amino acid residue."""
    aa_to_ccd = get_amino_acid_ccd_map()
    ligand_ccd = aa_to_ccd.get(residue_char, 'UNK')
    
    return {
        "ligand": {
            "id": f"{chain_id}{ligand_suffix}",
            "ccdCodes": [ligand_ccd]
        }
    }

def create_protein_sequence(chain_id: str, sequence: str) -> Dict:
    """Create a protein sequence entry."""
    return {
        "protein": {
            "id": chain_id,
            "sequence": sequence
        }
    }

def add_peptide_bond(new_bonded_pairs: List, chain1_id: str, chain1_pos: int, chain2_id: str, chain2_pos: int):
    """Add a peptide bond between two chains."""
    new_bonded_pairs.append([[chain1_id, chain1_pos, "C"], [chain2_id, chain2_pos, "N"]])

def process_intra_chain_bond(modified_json: Dict, bond: Tuple, original_json: Dict) -> Dict:
    """Process bonds within the same chain."""
    atom1, atom2 = bond
    chain_id_orig, seq_num1_old, atom_name1 = atom1
    _, seq_num2_old, atom_name2 = atom2
    
    # Get the current chain info
    chain_info = get_sequence_info(original_json, chain_id_orig)
    if not chain_info:
        print(f"Warning: Could not find sequence info for chain {chain_id_orig}")
        return original_json
    
    sequence = chain_info["sequence"]
    pos1, pos2 = sorted([seq_num1_old, seq_num2_old])
    
    # Adjust atom names if positions were swapped
    if seq_num1_old > seq_num2_old:
        atom_name1, atom_name2 = atom_name2, atom_name1
    
    new_sequences = []
    new_bonded_pairs = []
    
    if pos1 == pos2:
        # Same residue - split into before/ligand/after
        parts = [
            (sequence[:pos1-1], f"{chain_id_orig}A") if pos1 > 1 else None,
            (sequence[pos1-1], f"{chain_id_orig}L"),
            (sequence[pos1:], f"{chain_id_orig}B") if pos1 < len(sequence) else None
        ]
        
        prev_chain_id = None
        for i, part in enumerate(parts):
            if part is None:
                continue
            seq_part, chain_id = part
            
            if i == 1:  # Ligand
                new_sequences.append(create_ligand_from_residue(chain_id_orig, seq_part))
            else:  # Protein
                new_sequences.append(create_protein_sequence(chain_id, seq_part))
            
            # Add peptide bonds
            if prev_chain_id and i == 1:  # Bond to ligand
                add_peptide_bond(new_bonded_pairs, prev_chain_id, len(parts[0][0]), chain_id, 1)
            elif prev_chain_id and i == 2:  # Bond from ligand
                add_peptide_bond(new_bonded_pairs, prev_chain_id, 1, chain_id, 1)
            
            prev_chain_id = chain_id if seq_part else None
    
    else:
        # Different residues - create multiple ligands
        parts = [
            (sequence[:pos1-1], f"{chain_id_orig}A") if pos1 > 1 else None,
            (sequence[pos1-1], f"{chain_id_orig}L1"),
            (sequence[pos1:pos2-1], f"{chain_id_orig}B") if pos2 > pos1 + 1 else None,
            (sequence[pos2-1], f"{chain_id_orig}L2"),
            (sequence[pos2:], f"{chain_id_orig}C") if pos2 < len(sequence) else None
        ]
        
        prev_chain_id = None
        for i, part in enumerate(parts):
            if part is None:
                continue
            seq_part, chain_id = part
            
            if i in [1, 3]:  # Ligands
                suffix = "L1" if i == 1 else "L2"
                new_sequences.append(create_ligand_from_residue(chain_id_orig, seq_part, suffix))
            else:  # Protein
                new_sequences.append(create_protein_sequence(chain_id, seq_part))
            
            # Add peptide bonds between adjacent parts
            if prev_chain_id:
                if i in [1, 3]:  # To ligand
                    prev_len = len(parts[i-1][0]) if parts[i-1] else 1
                    add_peptide_bond(new_bonded_pairs, prev_chain_id, prev_len, chain_id, 1)
                else:  # From ligand
                    add_peptide_bond(new_bonded_pairs, prev_chain_id, 1, chain_id, 1)
            
            prev_chain_id = chain_id if seq_part else None
        
        # Add the original bond between ligands
        new_bonded_pairs.append([[f"{chain_id_orig}L1", 1, atom_name1], 
                               [f"{chain_id_orig}L2", 1, atom_name2]])
    
    # Add other sequences unchanged
    for sequence_item in original_json["sequences"]:
        if "protein" in sequence_item and sequence_item["protein"]["id"] != chain_id_orig:
            new_sequences.append(sequence_item)
        elif "ligand" in sequence_item or "dna" in sequence_item or "rna" in sequence_item:
            new_sequences.append(sequence_item)
    
    # Add existing bonds (excluding the original bond)
    for existing_bond in original_json.get("bondedAtomPairs", []):
        if tuple(existing_bond) != bond:
            new_bonded_pairs.append(existing_bond)
    
    modified_json["sequences"] = new_sequences
    modified_json["bondedAtomPairs"] = new_bonded_pairs
    return modified_json

def process_chain_bond(modified_json: Dict, bond: Tuple, is_intra_chain: bool,residue_mapping: Dict, original_json: Dict) -> Dict:
    """Process bonds between different chains."""
    atom1, atom2 = bond
    chain1_id_orig, seq_num1_old, atom_name1 = atom1
    chain2_id_orig, seq_num2_old, atom_name2 = atom2
    chain1_id = residue_mapping[chain1_id_orig][seq_num1_old]["modified_chain_id"]
    chain2_id = residue_mapping[chain2_id_orig][seq_num2_old]["modified_chain_id"]
    seq_num1 = residue_mapping[chain1_id_orig][seq_num1_old]["modified_residue_num"]
    seq_num2 = residue_mapping[chain2_id_orig][seq_num2_old]["modified_residue_num"]
    # Get sequence information
    chain1_info = get_sequence_info(modified_json, chain1_id)
    chain2_info = get_sequence_info(modified_json, chain2_id)
    
    if not chain1_info or not chain2_info:
        print(f"Warning: Could not find sequence info for chains {chain1_id_orig} or {chain2_id_orig}")
        return modified_json
    
    chain1_sequence = chain1_info["sequence"]
    chain2_sequence = chain2_info["sequence"]
    
    # Determine which residue to use as ligand (prefer terminal residues)
    chain1_is_terminal = is_terminal_residue(chain1_sequence, seq_num1)
    chain2_is_terminal = is_terminal_residue(chain2_sequence, seq_num2)
    use_chain1_as_ligand = chain1_is_terminal or not chain2_is_terminal
    
    new_sequences = []
    new_bonded_pairs = []
    
    if use_chain1_as_ligand:
        # Convert chain1 residue to ligand
        ligand_chain_id = chain1_id_orig
        ligand_seq_num = seq_num1_old
        ligand_residue = chain1_sequence[seq_num1_old - 1]
        ligand_atom_name = atom_name1
        target_chain_id = chain2_id_orig
        target_seq_num = seq_num2_old
        target_atom_name = atom_name2
        ligand_sequence = chain1_sequence
        is_ligand_terminal = chain1_is_terminal
    else:
        # Convert chain2 residue to ligand
        ligand_chain_id = chain2_id_orig
        ligand_seq_num = seq_num2_old
        ligand_residue = chain2_sequence[seq_num2_old - 1]
        ligand_atom_name = atom_name2
        target_chain_id = chain1_id_orig
        target_seq_num = seq_num1_old
        target_atom_name = atom_name1
        ligand_sequence = chain2_sequence
        is_ligand_terminal = chain2_is_terminal
    
    # Create ligand
    ligand_id = f"{ligand_chain_id}L"
    new_sequences.append(create_ligand_from_residue(ligand_chain_id, ligand_residue))
    
    # Handle ligand chain splitting
    if is_ligand_terminal:
        is_c_terminal = ligand_seq_num == len(ligand_sequence)
        update_residue_mapping_for_terminal_split(residue_mapping, ligand_chain_id, 
                                                ligand_seq_num, len(ligand_sequence), is_c_terminal)
        
        if is_c_terminal:
            modified_sequence = ligand_sequence[:-1]
            new_chain_id = f"{ligand_chain_id}A"
        else:
            modified_sequence = ligand_sequence[1:]
            new_chain_id = f"{ligand_chain_id}A"
        
        if modified_sequence:
            new_sequences.append(create_protein_sequence(new_chain_id, modified_sequence))
            
            # Add peptide bonds
            if is_c_terminal:
                add_peptide_bond(new_bonded_pairs, new_chain_id, len(modified_sequence), ligand_id, 1)
            else:
                add_peptide_bond(new_bonded_pairs, ligand_id, 1, new_chain_id, 1)
    else:
        # Internal residue - split into two parts
        update_residue_mapping_for_internal_split(residue_mapping, ligand_chain_id, 
                                                ligand_seq_num, len(ligand_sequence))
        
        part_a_sequence = ligand_sequence[:ligand_seq_num-1]
        part_b_sequence = ligand_sequence[ligand_seq_num:]
        
        if part_a_sequence:
            part_a_id = f"{ligand_chain_id}A"
            new_sequences.append(create_protein_sequence(part_a_id, part_a_sequence))
            add_peptide_bond(new_bonded_pairs, part_a_id, len(part_a_sequence), ligand_id, 1)
        
        if part_b_sequence:
            part_b_id = f"{ligand_chain_id}B"
            new_sequences.append(create_protein_sequence(part_b_id, part_b_sequence))
            add_peptide_bond(new_bonded_pairs, ligand_id, 1, part_b_id, 1)
    
    # Add other sequences unchanged (excluding the ligand chain)
    for sequence in modified_json["sequences"]:
        if "protein" in sequence and sequence["protein"]["id"] != ligand_chain_id:
            new_sequences.append(sequence)
        elif "ligand" in sequence or "dna" in sequence or "rna" in sequence:
            new_sequences.append(sequence)
    
    # Bond from ligand to target chain
    target_mapped = residue_mapping[target_chain_id][target_seq_num]
    new_bonded_pairs.append([[ligand_id, 1, ligand_atom_name], 
                           [target_mapped["modified_chain_id"], 
                            target_mapped["modified_residue_num"], target_atom_name]])
    
    # Add existing bonds (excluding the original bond)
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
        print(f"\nProcessing {filename}...")
        
        # Initialize residue mapping
        residue_mapping = initialize_residue_mapping(json_data)
        
        # Find protein-protein bonds
        protein_bonds = find_protein_protein_bonds(json_data)
        
        if not protein_bonds:
            print(f"No protein-protein bonds found in {filename}")
            continue
        
        # Process each protein-protein bond
        modified_json = json_data
        for i, bond in enumerate(protein_bonds, 1):
            print(f"Modifying bond {i}/{len(protein_bonds)}: {bond}")
            modified_json = model_bond_with_ligand(modified_json, bond, residue_mapping)
        
        # Save modified JSON
        output_path = Path(output_dir) / f"{filename}_modified.json"
        with open(output_path, 'w') as f:
            json.dump(modified_json, f, indent=2)
        
        print(f"Saved modified file: {output_path}")

def main():
    """
    Main function to execute the protein bond modeling script.
    """
    parser = argparse.ArgumentParser(
        description="Model protein-protein bonds using ligand bridges in AlphaFold3 JSON files"
    )
    parser.add_argument(
        "--source-dir", 
        "-s", 
        default="test_files/input/",
        help="Directory containing input JSON files (default: output/jsons/ubn_links/)"
    )
    parser.add_argument(
        "--output-dir", 
        "-o", 
        default="test_files/output/",
        help="Directory to save modified JSON files (default: output/jsons/ubn_links_modified/)"
    )
    parser.add_argument(
        "--verbose", 
        "-v", 
        action="store_true",
        help="Enable verbose output"
    )
    
    args = parser.parse_args()
    
    print("Starting protein bond modeling process...")
    print(f"Source directory: {args.source_dir}")
    print(f"Output directory: {args.output_dir}")
    
    if args.verbose:
        print("Verbose mode enabled")
    
    # Check if source directory exists
    if not Path(args.source_dir).exists():
        print(f"Error: Source directory '{args.source_dir}' does not exist!")
        return 1
    
    # Process all JSON files
    process_json_files(args.source_dir, args.output_dir)
    print("Process completed successfully!")
    return 0

if __name__ == "__main__":
    exit_code = main()
    exit(exit_code)