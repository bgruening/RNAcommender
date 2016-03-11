#!/usr/bin/env python

from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import sys
from eden.converter.fasta import fasta_to_sequence
from eden.converter.rna.rnaplfold import rnaplfold_to_eden
from eden.graph import Vectorizer
from eden.util import vectorize as eden_vectorize
from sklearn import preprocessing

import fasta_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

class RNAVectorizer():
    """Compute the RNA features"""
    def __init__(self,fasta,output,window_size=150,max_bp_span=40,avg_bp_prob_cutoff=0.4,
        complexity=2,nbits=10,njobs=-1,verbose=True):
        """
        Params
        ------
        fasta : str
            Fasta file containing the RNA sequences.

        output : str
            Name of the output file. The output file is an HDF5 containing a
            pandas DataFrame, in which the columns are the RNA names and the rows
            are the EDeN features.

        window_size : int (default : 150)
            Window size of RNAplfold. Average the pair
            probabilities over windows of given size.

        max_bp_span : int (default : 40)
            Maximum allowed separation of a base pair to span.
            I.e. no pairs (i,j) with j-i > span will be allowed.

        avg_bp_prob_cutoff : float (default : 0.4)
            Report only base pairs with an average probability > cutoff.

        complexity : int (default : 2)
            Complexity of the features extracted. Equivalent to
            define EDeN parameters d = r = complexity.

        nbits : int (default : 10)
            Number of bits that defines the feature space size: |feature space|=2^nbits.

        njobs : int (default : -1)
            Number of parallel jobs (default: all CPUs).

        verbose : bool (default : True)
            Print information to STDOUT.
        """
        self.fasta = fasta
        self.output = output
        self.window_size=window_size
        self.max_bp_span = max_bp_span
        self.avg_bp_prob_cutoff=avg_bp_prob_cutoff
        self.complexity = complexity
        self.nbits = nbits
        self.njobs = njobs
        self.verbose = verbose

    def _fold_sequences(self):
        """ Fold the RNA sequences using RNAplfold"""
        if self.verbose:
            print("Folding sequences using RNAplfold -W %i -L %i -c %f --noLP..." % (self.window_size,self.max_bp_span,self.avg_bp_prob_cutoff), end= ' ')
            sys.stdout.flush()

        seqs = fasta_to_sequence(self.fasta)
        graphs = rnaplfold_to_eden(seqs,
                               window_size=self.window_size,
                               max_bp_span=self.max_bp_span,
                               avg_bp_prob_cutoff=self.avg_bp_prob_cutoff,
                               max_num_edges=1)
        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()
        return graphs

    def _vectorize_graphs(self,graphs):
        """Vectorize the RNAplfold graphs using EDeN"""
        if self.verbose:
            print("Vectorizing (complexity: %i, hashing: %i bits)..." % (self.complexity,self.nbits), end= ' ')
            sys.stdout.flush()

        vec = Vectorizer(complexity=self.complexity,nbits=self.nbits)
        X_sparse = eden_vectorize(graphs, vectorizer=vec, n_jobs=self.njobs)

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()
        return X_sparse.todense()

    def vectorize(self):
        """Produce the RNAfeatures"""
        names = fasta_utils.seq_names(self.fasta)

        graphs = self._fold_sequences()
        X = self._vectorize_graphs(graphs)

        df = pd.DataFrame(X.T[1:],columns=names)
        store = pd.io.pytables.HDFStore(self.output)
        store['features'] = df
        store.close()

        if self.verbose:
            print("RNA features saved in %s" % self.output)
            sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta', metavar='fasta', type=str,
                        help="""Fasta file containing the RNA sequences.""")
    parser.add_argument('output', metavar='output', type=str,
                        help="""File name of the HDF Store to save the RNA features.""")
    # RNAplfold parameters
    parser.add_argument('--window-size', metavar='window_size', type=int, default=150,
                        help="""Window size of RNAplfold.""")
    parser.add_argument('--max-bp-span', metavar='max_bp_span', type=int, default=40,
                        help="""Maximum allowed separation of a base pair to span.""")
    parser.add_argument('--avg-bp-prob-cutoff', metavar='avg_bp_prob_cutoff', type=float, default=0.4,
                        help="""Report only base pairs with an average probability > cutoff.""")
    # EDeN parameters
    parser.add_argument('--complexity', metavar='complexity', type=int, default=2,
                        help="""Complexity of the features extracted.""")
    parser.add_argument('--nbits', metavar='nbits', type=int, default=10,
                        help="""Number of bits that defines the feature space size: |feature space|=2^nbits.""")
    # Other paramentes
    parser.add_argument('--njobs', metavar='njobs', type=int, default=-1,
                        help="""Number of parallel jobs (-1 means all CPUs).""")
    parser.add_argument('--quiet', dest='quiet', action='store_true', default=False,
                        help="""Do not print information at STDOUT.""")

    args = parser.parse_args()

    v = RNAVectorizer(fasta=args.fasta,
                    output=args.output,
                    window_size=args.window_size,
                    max_bp_span=args.max_bp_span,
                    avg_bp_prob_cutoff=args.avg_bp_prob_cutoff,
                    complexity=args.complexity,
                    nbits=args.nbits,
                    njobs=args.njobs,
                    verbose=(not args.quiet))
    v.vectorize()
