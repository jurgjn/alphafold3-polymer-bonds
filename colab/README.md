## Run AlphaFold3 on Google Colab (Pro)
Install AlphaFold3 inference pipeline on Google Colab
- [Built-in unit tests pass](https://colab.research.google.com/github/jurgjn/alphafold3-polymer-bonds/blob/master/colab/run_alphafold3_tests.ipynb)
- [Run full inference, e.g any ligand](https://colab.research.google.com/github/jurgjn/alphafold3-polymer-bonds/blob/master/colab/run_alphafold3_inference.ipynb)

Known caveats:
- Not enough disk space for the sequence databases, need to use MSAs from elsewhere
- Limited input size, e.g. ~3k tokens on an A100 with 40 GB GPU RAM
