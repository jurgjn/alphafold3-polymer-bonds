# NOTES
- Coordinates in AlphaFold3 input are 1-based (e.g. `ptmType`, `bondedAtomPairs`)
- To introduce an atom-level representation of a residue for covalent bond definition, use `ptmType`/`userCcd` to define a ligand that's identical to the standard amino acid
[alphafold3:159](https://github.com/google-deepmind/alphafold3/issues/159#issuecomment-2525070335)
- Leaving atoms for polymer-ligand covalent bonds were kept when training AlphaFold3 
[alphafold3:159](https://github.com/google-deepmind/alphafold3/issues/159#issuecomment-2523711478)
- Cif is heterogenous format, mandatory fields for AlphaFold3 ligands are documented in
[alphafold3:#178](https://github.com/google-deepmind/alphafold3/issues/178#issuecomment-2521175288)
- Modified residues need to have the following atoms: `N`, `CA`, `C`, `O`, `OXT`
[alphafold3:#159](https://github.com/google-deepmind/alphafold3/issues/159#issuecomment-2561311898)
- Cif files can include information on leaving atoms but AlphaFold3 ignores this
[alphafold3:#250](https://github.com/google-deepmind/alphafold3/issues/159#issue-2712293489)
- Leaving atoms are harmless as the model does not seem to use them
[alphafold3:#250](https://github.com/google-deepmind/alphafold3/issues/250#issuecomment-2580322870)
- Warnings about `does not contain a pseudo-beta atom` described as safe to ignore
[alphafold3:#438](https://github.com/google-deepmind/alphafold3/issues/438#issuecomment-2955474005)
- AlphaFold server allows for modifications for
[residues](https://github.com/google-deepmind/alphafold/tree/main/server#protein-chains),
[DNA](https://github.com/google-deepmind/alphafold/tree/main/server#dna-chains),
[RNA](https://github.com/google-deepmind/alphafold/tree/main/server#rna-chains)
