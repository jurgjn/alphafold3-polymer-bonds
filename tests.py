import json
import unittest
import subprocess
import sys
from pathlib import Path

class TestBondModeling(unittest.TestCase):
    
    def setUp(self):
        """Set up test directories."""
        self.input_dir = Path("test_files/input")
        self.output_dir = Path("test_files/output")
        self.solution_dir = Path("test_files/solution")
        
        # Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def compare_json_files(self, file1_path, file2_path):
        """Compare two JSON files ignoring dictionary order."""
        try:
            with open(file1_path, 'r') as f1:
                json1 = json.load(f1)
            with open(file2_path, 'r') as f2:
                json2 = json.load(f2)
            
            # Compare the JSON objects (order-independent)
            return json1 == json2
        except Exception as e:
            print(f"Error comparing {file1_path} and {file2_path}: {e}")
            return False
    
    def test_json_files(self):
        """Test all JSON files in input folder against solution folder."""
        # First, run the protein bond modeling script
        print("Running protein bond modeling script...")
        
        try:
            result = subprocess.run([
                sys.executable, "model_protein_bonds_hack.py",
                "--source-dir", "test_files/input",
                "--output-dir", "test_files/output",
                "--verbose"
            ], capture_output=True, text=True, check=True)
            
            print("Script output:")
            print(result.stdout)
            if result.stderr:
                print("Script errors:")
                print(result.stderr)
                
        except subprocess.CalledProcessError as e:
            self.fail(f"Protein bond modeling script failed: {e}\nStdout: {e.stdout}\nStderr: {e.stderr}")
        
        # Get all input files
        input_files = list(self.input_dir.glob("*.json"))
        self.assertGreater(len(input_files), 0, "No input JSON files found")
        
        print(f"\nTesting {len(input_files)} files...")
        
        # Check each processed file against its solution
        for input_file in input_files:
            with self.subTest(file=input_file.name):
                # Expected output file name (with _modified suffix)
                output_file = self.output_dir / f"{input_file.stem}_modified.json"
                solution_file = self.solution_dir / f"{input_file.stem}_modified.json"
                
                # Check if output file was created
                self.assertTrue(output_file.exists(), 
                              f"Output file {output_file} was not created")
                
                # Check if solution file exists
                self.assertTrue(solution_file.exists(), 
                              f"Solution file {solution_file} does not exist")
                
                # Compare the files
                files_match = self.compare_json_files(output_file, solution_file)
                
                self.assertTrue(files_match, 
                              f"Generated file {output_file.name} does not match solution {solution_file.name}")
        
        print(f"\nAll {len(input_files)} test files passed!")

if __name__ == "__main__":
    unittest.main()