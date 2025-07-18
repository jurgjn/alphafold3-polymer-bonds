# AF3 Protein-Protein Bond Modeling Tool

A Python tool for modeling protein-protein bonds in AlphaFold3 structures using ligand bridges. This tool processes JSON files containing protein structure data and converts direct protein-protein bonds into chemically realistic ligand-mediated connections.

## Overview

This tool addresses the challenge of representing protein-protein bonds in structural biology by:
- Converting amino acids involved in protein-protein bonds into ligand molecules
- Splitting protein chains appropriately to maintain proper connectivity
- Creating chemically realistic bond networks through ligand intermediates

## Quick Start

### Prerequisites
- Python 3.8+
- Input JSON files with AlphaFold3 structure data

### Installation
```bash
git clone <repository-url>
cd model-protein-bond
```

### Basic Usage
```bash
# Process all JSON files in the input directory
python model_protein_bonds_hack.py --source-dir test_files/input/ --output-dir test_files/output/

# With verbose output
python model_protein_bonds_hack.py -s test_files/input/ -o test_files/output/ --verbose
```

### File Structure
```
model-protein-bond/
├── model_protein_bonds_hack.py    # Main processing script
├── README.md                      # This file
├── tests.py                       # Test suite
└── test_files/
    ├── input/                     # Input JSON files
    ├── output/                    # Generated output files
    └── solution/                  # Reference solutions
```

## Command Line Options

| Option | Short | Description | Default |
|--------|-------|-------------|---------|
| `--source-dir` | `-s` | Directory containing input JSON files | `test_files/input/` |
| `--output-dir` | `-o` | Directory for output files | `test_files/output/` |
| `--verbose` | `-v` | Enable detailed output | `False` |

## Example Usage
```bash
# Process files from custom directories
python model_protein_bonds_hack.py --source-dir input/ --output-dir output/
python model_protein_bonds_hack.py -s input/ -o output/ --verbose
```

## How It Works

### Algorithm Overview
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

## Input/Output Format

### Input Files
- **Format**: JSON files containing AlphaFold3 structure data
- **Required fields**: 
  - `sequences`: Array of protein/ligand sequence definitions
  - `bondedAtomPairs`: Array of atomic bonds in the structure
- **Location**: Specified by `--source-dir` parameter

### Output Files
- **Format**: Modified JSON files with ligand bridge representations
- **Naming**: Original filename with `_modified.json` suffix
- **Location**: Specified by `--output-dir` parameter

### Example Transformation
```
Original: protein_structure.json
Output:   protein_structure_modified.json
```

## Chain ID Conventions

When the tool processes protein-protein bonds, it creates new chain identifiers:

| Original | After Processing | Description |
|----------|------------------|-------------|
| `A` | `AA`, `AL`, `AB` | Chain A split with ligand L between parts A and B |
| `B` | `BA`, `BL` | Chain B terminal residue becomes ligand L |

## Testing

Run the test suite to verify functionality:
```bash
python tests.py
```

Tests compare generated output files with reference solutions in `test_files/solution/`.