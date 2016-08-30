# RNAcommender
RNAcommender is a tool for genome-wide recommendation of RNA-protein interactions. It is a recommender system capable of suggesting RNA targets to unexplored RNA binding proteins, by propagating the available interaction information, taking into account the protein domain composition and the RNA predicted secondary structure.


**G. Corrado, T. Tebaldi, F. Costa, P. Frasconi and A. Passerini. "RNAcommender: genome-wide recommendation of RNA-protein interactions." *Bioinformatics*, 2016.**

Requirements
------------
Include the following executable files in your PATH:
* From SAM (3.5) [[download](https://compbio.soe.ucsc.edu/sam2src/)]:
    - build_model
    - get_fisher_scores
* From the Vienna RNA package (2.1.8) [[download](https://www.tbi.univie.ac.at/RNA/)]:
    - RNAplfold

* Python packages required (and versions used during development):
    - numpy 1.8.1
    - pandas 0.17.1
    - tables 3.2.2.dev0
    - Theano 0.7

* The RNAfeatures are generated using EDeN. You can install EDeN with pip directly from github.

```bash
pip install git+https://github.com/fabriziocosta/EDeN.git --user
```

Training
========
Here the instructions to train a model using RNAcommender. Complete documentation of each python script can be accessed using the -h command. In this package we provide a recommender trained, using cross-validated parameters, from the high-thoughput human interaction in the AURA 2.5 dataset (```example/trained_from_AURA_HT.pkl```). This trained model can be used to predict protein-RNA interactions. If you want to use this trained model to get protein-RNA interactions recommendations you can skip the Training section, and move to the Recommending one.

Protein features (training)
---------------------------
```rbpfeatures.py``` produces the protein features. It requires two fasta files: one with the sequences of the proteins used as reference for the similarity, and one for the proteins for which we want to have the features. When preparing the protein features for the training we want these two set to be exactly identical and to contain only the proteins present in the training set. We compute the features for the training proteins by executing the following command:

```bash
python rbpfeatures.py ../examples/rbps_HT.fa ../examples/rbps_HT.fa ../examples/rbps_HT.h5 --all-sel
```
rbps_HT.h5 is the output file that will store the features, the flag ```--all-sel``` forces the inclusion of all the selected sequences in the output. In this case if a train protein has no similarity with the other proteins we still want it to use it (this protein will be represented by its one-hot encoding).

RNA features
------------
```rnafeatures.py``` produces the RNA features. It requires in input the fasta file and the name of the output file that will store the features:

```bash
python rnafeatures.py ../examples/utrs.fa  ../examples/utrs.h5
```

Interaction matrix
------------------
```interaction.py``` prepares the interaction matrix. It requires a file containing the interaction map (see interactions_HT.txt in examples), and the name of the output file where to store the interaction matrix:

```bash
python interactions.py ../examples/interactions_HT.txt  ../examples/interactions_HT.h5
```

Training the recommender
------------------------
At this point we have all the information required to train a model. ```train.py``` allows to train a recommender that later will be used to predict protein-RNA interactions. In order to train the model we need to specify the protein features, the RNA features, the interaction matrix, and where to save the trained model:

```bash
python train.py ../examples/rbps_HT.h5 ../examples/utrs.h5 ../examples/interactions_HT.h5 ../examples/trained_recommender.pkl --standardize-Fr
```
For the origin of the features we want to standardize the RNA features only, and we do that by activating the flag ```--standardize-Fr```.
Run ```python train.py -h``` for the list of parameters of the model. Parameter settings vary from dataset to dataset and we encourage to use cross-validated parameters.

Recommending
============
Here the instructions to train a model using RNAcommender. Complete documentation of each python script can be accessed using the -h command. RNAcommender produces a ranked list of protein-RNA interactions.

NOTE: if you want to use our pretrained recommender be sure to set the Theano flag floatX=float32.

Protein features (recommending)
-------------------------------
```rbpfeatures.py``` produces the protein features. It requires two fasta files: one with the sequences of the proteins used as reference for the similarity (which are the one used for training the model), and one for the proteins for which we want to have the features (which are the unexplored proteins). We compute the features for the unexplored proteins by executing the following command:

```bash
python rbpfeatures.py ../examples/rbps_HT.fa ../examples/rbps_new.fa ../examples/rbps_new.h5
```
rbps_new.h5 is the output file that will store the features. In this case the flag ```--all-sel``` MUST NOT be used. We need to discard all the unknown proteins that have zero similarity with the proteins in the training set, because, for them, it is not possible to perform de novo recommendations.

RNA features
------------
RNA features are computed in the same way as explained in the Training section. ```rnafeatures.py``` produces the RNA features. It requires in input the fasta file and the name of the output file that will store the features:

```bash
python rnafeatures.py ../examples/utrs.fa  ../examples/utrs.h5
```

Obtaining the recommendations
-----------------------------
```recommend.py``` uses a trained recommender to return a ranked list of protein-RNA interactions. It requires the features for the unknown proteins, and the features for the RNAs and the trained model. Additionally we can specify an output file to store the results (otherwise they will be print at STDOUT):

```bash
python recommend.py ../examples/rbps_new.h5  ../examples/utrs.h5 ../examples/trained_from_AURA_HT.pkl --output ../examples/recommendations.txt --standardize-Fr
```
We added again the flag ```--standardize-Fr``` to match the case used during training.

It is also possible to specify one or more protein that will be included in the results (all the others will be discarded). For example:
```bash
python recommend.py ../examples/rbps_new.h5  ../examples/utrs.h5 ../examples/trained_from_AURA_HT.pkl --output ../examples/recommendations.txt --standardize-Fr --to-predict RALY
```
will recommend RNA targets only to the protein RALY, while:
```bash
python recommend.py ../examples/rbps_new.h5  ../examples/utrs.h5 ../examples/trained_from_AURA_HT.pkl --output ../examples/recommendations.txt --standardize-Fr --to-predict ANKHD1 RALY
```
will recommend RNA targets to the proteins ANKHD1 and RALY.
