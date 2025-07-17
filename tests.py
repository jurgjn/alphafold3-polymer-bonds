import json
import unittest
from pathlib import Path

def model_bonded_atoms(input_data):
    """Main function that processes bonded atoms."""
    # Import your functions from the notebook
    from your_module import initialize_residue_mapping, find_protein_protein_bonds, model_bond_with_ligand
    
    residue_mapping = initialize_residue_mapping(input_data)
    bonds = find_protein_protein_bonds(input_data)
    
    result = input_data
    for bond in bonds:
        result = model_bond_with_ligand(result, bond, residue_mapping)
    
    return result

class TestBondModeling(unittest.TestCase):
    
    def setUp(self):
        """Set up test directories."""
        self.input_dir = Path("test_files/input")
        self.solution_dir = Path("test_files/solution")
    
    def test_json_files(self):
        """Test all JSON files in input folder against solution folder."""
        # Get all input files
        input_files = list(self.input_dir.glob("*.json"))
        
        for input_file in input_files:
            with self.subTest(file=input_file.name):
                # Load input
                with open(input_file, 'r') as f:
                    input_data = json.load(f)
                
                # Load expected solution
                solution_file = self.solution_dir / input_file.name
                with open(solution_file, 'r') as f:
                    expected = json.load(f)
                
                # Run the function
                result = model_bonded_atoms(input_data)
                
                # Compare
                self.assertEqual(result, expected)

if __name__ == "__main__":
    unittest.main()