[![PyPI - Version](https://img.shields.io/pypi/v/alphafold3-polymer-bonds)](https://pypi.org/project/alphafold3-polymer-bonds/)
[![PyPI Downloads](https://static.pepy.tech/personalized-badge/alphafold3-polymer-bonds?period=total&units=INTERNATIONAL_SYSTEM&left_color=BLACK&right_color=GREEN&left_text=downloads)](https://pepy.tech/projects/alphafold3-polymer-bonds)

# Polymer bonds in AlphaFold3
AlphaFold3
[does not allow](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md#bonds)
covalent bonds between/within polymer chains (protein, DNA, RNA).
We work around this limitation by treating one of the corresponding residues/nucleic acids as a ligand.
In principle, this may enable AlphaFold3 to explicitly model e.g. disulfide bonds, cyclic peptides, zero-length crosslinkers, protein-DNA bonds..

*This is currently exploratory, see below for specifc examples. Also take a look at complementary work:
[KosinskiLab/af3x](https://github.com/KosinskiLab/af3x)
and
[bio-phys/polyUb-AF](https://github.com/bio-phys/polyUb-AF).*

![1DF6](https://raw.githubusercontent.com/jurgjn/alphafold3-polymer-bonds/main/examples/visualise/1DF6.png)

*[1DF6](https://www.rcsb.org/structure/1DF6): macrocyclic peptide cycloviolacin O1 with (left, RMSD=1.979) and without (right, RMSD=2.632) a covalent bond between Ser1 and Glu30. The model with the explicit covalent bond reproduces three cysteine bridges as in the experimental structure. The model without the covalent bond does not.*

![6OQ1](https://raw.githubusercontent.com/jurgjn/alphafold3-polymer-bonds/main/examples/visualise/6OQ1.png)

*[6OQ1](https://www.rcsb.org/structure/6OQ1): K11/K48 branched tri-ubiquitin with (left, RMSD=6.157) and without (right, RMSD=9.676) covalent bonds between K11/K48 (orange side chains) and Gly76 (red side chains).*

![8S6W](https://raw.githubusercontent.com/jurgjn/alphafold3-polymer-bonds/main/examples/visualise/8S6W.png)

*[8S6W](https://www.rcsb.org/structure/8S6W): RNA complex with two circular chains with (left, RMSD=9.788) and without (right, RMSD=23.233) covalent bonds to define the circularisation. The model with the explicit circularisation coarsely reproduces the overall shape of the experimental structure. The model without circularisation bonds does not.*

## Quick start
Install using `pip`, the main script is called `alphafold3_polymer_bonds`. Use `--bonds_path` to specify a standard AlphaFold 3 input .json with extra polymer bonds in `bondedAtomPairs`. The script will output a modified .json (`--encoded_path encoded_input.json`) that can be used as input for AlphaFold 3.

```bash
pip install alphafold3-polymer-bonds
alphafold3_polymer_bonds --bonds_path input_with_polymer_bonds.json --encoded_path encoded_input.json
```
