# RNAcommender
RNAcommender is a tool for genome-wide recommendation of RNA-protein interactions. It is a recommender system capable of suggesting RNA targets to unexplored RNA binding proteins, by propagating the available interaction information, taking into account the protein domain composition and the RNA predicted secondary structure.

Requirements
------------
Include the following executable files in your PATH:
* From SAM (3.5) [[download](https://compbio.soe.ucsc.edu/sam2src/)]:
    - build_model
    - get_fisher_scores
* From the Vienna RNA package (2.1.8) [[download](https://www.tbi.univie.ac.at/RNA/)]:
    - RNAplfold

* Python packages required (and versions used during test):
    - numpy 1.8.1
    - pandas 0.17.1
    - tables 3.2.2.dev0
    - scikit-learn 0.18.dev0
    - Theano 0.7

* The RNAfeatures are generated using EDeN. You can install EDeN with pip directly from github.

```python
pip install git+https://github.com/fabriziocosta/EDeN.git --user
```

