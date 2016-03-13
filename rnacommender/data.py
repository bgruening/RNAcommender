from __future__ import print_function
import sys
import numpy as np
import pandas as pd
from theano import config

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

class Dataset(object):
    def __init__(self,Fp,Fr,
        standardize_proteins=False,standardize_rnas=False,
        verbose=True):
        """
        Parameters
        ----------
        Fp : str
            The name of the HDF5 file containing features for the proteins.

        Fr : str
            The name of the HDF5 file containing features for the RNAs.

        standardize_proteins : bool (default : False)
            Whether protein features should be standardized.

        standardize_rnas : bool (default : False)
            Whether RNAs features should be standardized.

        verbose : bool (default : True)
            Print information at STDOUT.
        """
        self.verbose = verbose

        def standardize(X):
            Xmean = X.mean(axis=1)
            Xstd = X.std(axis=1)
            Xstd[Xstd == 0.0] = 1.0
            return ((X.T - Xmean) / Xstd).T

        store = pd.io.pytables.HDFStore(Fp)
        self.Fp = store.features.astype(config.floatX)
        store.close()
        if self.verbose:
            print('Protein features of shape', self.Fp.shape)
            sys.stdout.flush()

        if standardize_proteins:
            if self.verbose:
                print('Standardizing protein features...', end = ' ')
                sys.stdout.flush()
            self.Fp = standardize(self.Fp)
            if self.verbose:
                print('Done.')
                sys.stdout.flush()

        store = pd.io.pytables.HDFStore(Fr)
        self.Fr = store.features.astype(config.floatX)
        store.close()
        if self.verbose:
            print('RNA features of shape', self.Fr.shape)
            sys.stdout.flush()

        if standardize_rnas:
            if self.verbose:
                print('Standardizing RNA features...', end = ' ')
                sys.stdout.flush()
            self.Fr = standardize(self.Fr)
            if self.verbose:
                print('.Done.')
                sys.stdout.flush()

    def load(self):
        raise NotImplementedError()

class TrainDataset(Dataset):
    def __init__(self,Fp,Fr,Y,
        standardize_proteins=False,standardize_rnas=False,
        verbose=True, seed=1234):
        """
        Parameters
        ----------
        Fp : str
            The name of the HDF5 file containing features for the proteins.

        Fr : str
            The name of the HDF5 file containing features for the RNAs.

        Y : str
            The name of the HDF5 file containing the interaction matrix.
            The dataframe inside this file has proteins as columns and RNAs as rows (index).
            A NaN is expected there if the interaction is supposedly unknown.

        standardize_proteins : bool (default : False)
            Whether protein features should be standardized.

        standardize_rnas : bool (default : False)
            Whether RNAs features should be standardized.

        verbose : bool (default : True)
            Print information at STDOUT.

        seed : int (default : 1234)
            Seed for random number generator.
        """
        super(TrainDataset,self).__init__(Fp,Fr,standardize_proteins,standardize_rnas,verbose)

        store = pd.io.pytables.HDFStore(Y)
        self.Y = store.matrix.astype(config.floatX)
        store.close()
        if self.verbose:
            print('Interaction matrix of shape', self.Y.shape)
            sys.stdout.flush()

        assert self.Fp.shape[1] == self.Y.shape[1] and self.Fr.shape[1] == self.Y.shape[0]

        self.seed = seed


    def load(self):
        """
        Returns
        -------
        dataset : list
            List of triplets (P,R,I) representing the batches.
            Each batch is made of all the labeled examples of one RNA.
        """
        if self.verbose:
            print('\nmaking training set (%d user%s and %d item%s)...' % (len(self.Y.columns), (len(self.Y.columns)>1)*'s',len(self.Y.index),(len(self.Y.index)>1)*'s'))
            sys.stdout.flush()

        protein_input_dim = self.Fp.shape[0]
        rna_input_dim = self.Fr.shape[0]
        dataset = []

        num_pos = 0
        num_neg = 0
        num_batches = 0

        progress = 0
        for (i,rna) in enumerate(self.Y.index):
            if i % (len(self.Y.index)/10) == 0:
                if self.verbose:
                    print(str(progress) + "%",end=' ')
                    sys.stdout.flush()
                progress += 10
            num_examples = self.Y.loc[rna].count().sum()  # understands NaN
            P = np.zeros((num_examples,protein_input_dim)).astype(config.floatX)
            R = np.zeros((num_examples,rna_input_dim)).astype(config.floatX)
            I = np.zeros((num_examples,1)).astype(config.floatX)

            index = 0
            for protein in self.Y.columns:
                if np.isnan(self.Y[protein][rna]):
                    continue
                P[index] = self.Fp[protein]
                R[index] = self.Fr[rna]
                I[index] = self.Y[protein][rna]
                if I[index] > 0:
                    num_pos += 1
                else:
                    num_neg += 1
                index += 1

            perm = np.random.permutation(range(num_examples))
            P = np.matrix(P[perm])
            R = np.matrix(R[perm])
            I = np.array(I[perm].flatten())

            dataset.append((P,R,I))
            num_batches += 1

        np.random.seed(self.seed)
        np.random.shuffle(dataset)

        if self.verbose:
            print("")
            print("Training set created")
            print("\twith %i examples" % (num_pos+num_neg))
            print("\tnumber of positives: %i" % num_pos)
            print("\tnumber of negatives: %i" % num_neg)
            print("\tnumber of batches: %i" % num_batches)
            sys.stdout.flush()

        return dataset

class PredictDataset(Dataset):
    def __init__(self,Fp,Fr,to_predict=None,
        standardize_proteins=False,standardize_rnas=False,
        verbose=True):
        """
        Parameters
        ----------
        Fp : str
            The name of the HDF5 file containing features for the proteins.

        Fr : str
            The name of the HDF5 file containing features for the RNAs.

        to_predict : list (default : None)
            List of proteins from Fp to predict.
            If None all the proteins will be predicted.

        standardize_proteins : bool (default : False)
            Whether protein features should be standardized.

        standardize_rnas : bool (default : False)
            Whether RNAs features should be standardized.

        verbose : bool (default : True)
            Print information at STDOUT.
        """
        super(PredictDataset,self).__init__(Fp,Fr,standardize_proteins,standardize_rnas,verbose)
        self.to_predict = to_predict
        self.Fp = self.Fp[self.to_predict]

    def load(self):
        """
        Returns
        -------
        Examples to predict. For each example:
            - P contains the protein features,
            - R contains the RNA features,
            - P_names contains the name of the protein,
            - R_names contains the name of the RNA.

        """
        if self.verbose:
            print('\nPreparing dataset (%d protein%s and %d RNA%s)...' % (self.Fp.shape[1], (self.Fp.shape[1]>1)*'s',self.Fr.shape[1],(self.Fr.shape[1]>1)*'s'), end=' ')
            sys.stdout.flush()

        protein_input_dim = self.Fp.shape[0]
        rna_input_dim = self.Fr.shape[0]
        num_examples = self.Fp.shape[1]*self.Fr.shape[1]
        P = np.zeros((num_examples,protein_input_dim)).astype(config.floatX)
        P_names = []
        R = np.zeros((num_examples,rna_input_dim)).astype(config.floatX)
        R_names = []
        index = 0
        for protein in self.Fp.columns:
            for rna in self.Fr.columns:
                P[index] = self.Fp[protein]
                P_names.append(protein)
                R[index] = self.Fr[rna]
                R_names.append(rna)
                index += 1

        if self.verbose:
            print("Done.\n")

        return(P,np.array(P_names),R,np.array(R_names))
