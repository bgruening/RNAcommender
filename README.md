# RNAcommender
RNAcommender is a tool for genome-wide recommendation of RNA-protein interactions. It is a recommender system capable of suggesting RNA targets to unexplored RNA binding proteins, by propagating the available interaction information, taking into account the protein domain composition and the RNA predicted secondary structure.

Requirements
============
RNAcommender requires some softwares that are not available as Python packages:
* From the Vienna RNA package (v 2.1.8 or later):
    - RNAplfold
* From SAM 3.5:
    - build_model
    - get_fisher_scores
