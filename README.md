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

```bash
pip install git+https://github.com/fabriziocosta/EDeN.git --user
```

Example
=======
Here the instructions to use RNAcommender on the datasets provided as examples. Complete documentation of each python script can be accessed using the -h command.

Protein features (training)
---------------------------
'''rbpfeatures.py''' produces the protein features. It requires two fasta files: one with the sequences of the proteins used as reference for the similarity, and one for the proteins for which we want to have the features. When computing the protein features for the protein in the training set we want these two set to be exactly identical. We compute the features for the training proteins by executing the following command:

```bash
python rbpfeatures.py ../examples/rbps_HT.fa ../examples/rbps_HT.fa ../examples/rbps_HT.h5 --all-sel
```
rbps_HT.h5 is the output file that will store the features, the flag --all-sel forces the inclusion of all the selected sequences in the output. In this case if a train protein has no similarity with the other proteins we still want it to use it (this protein will be represented by its one-hot encoding).

RNA features
------------

