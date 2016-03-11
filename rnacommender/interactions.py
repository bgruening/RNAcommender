#!/usr/bin/env python

from __future__ import print_function
import sys
import argparse
import numpy as np
import pandas as pd

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

class InteractionMatrix():
    """Prepare interaction matrix from interaction list"""
    def __init__(self,interaction_list,output,low_throughput_rbps=None,verbose=True):
        """
        Params
        ------
        interaction_list : str
            File containing the interaction list. Comments and header
            must start with #.

        output : str
            Name of the output file. The output file is an HDF5 containing a
            pandas DataFrame inside this file has proteins as columns and RNAs as rows (index).

        low_throughput_rbps : str (default : None)
            File containing a list (new line separated) of RBPs
            with low-throughput evidence only. Their missing interactions
            will be labeled with NaN.

        verbose : bool (default : True)
            Print information to STDOUT.
        """
        self.interaction_list = interaction_list
        self.output = output
        if low_throughput_rbps is not None:
            f = open(low_throughput_rbps)
            self.low_throughput_rbps = f.read().strip().split()
            f.close()
        else:
            self.low_throughput_rbps = []
        self.verbose = verbose

    def prepare(self):
        if self.verbose:
            print("Preparing interaction matrix for", end=' ')
        #read interaction list
        data = pd.read_table(self.interaction_list,sep='\t',comment='#',header=None)
        # determine the names of the RBPs and the RNAs
        rbps = data[0].unique()
        rnas = data[1].unique()
        # Define the DataFrame with all the intarction at 0
        n = rbps.shape[0]
        m = rnas.shape[0]

        if self.verbose:
            print("%i proteins and %i RNAs..." % (n,m), end=' ')

        Y = pd.DataFrame(np.zeros((m,n)),columns=rbps,index=rnas)
        # Set NaN the interactions for the low-throughput RBPs
        for rbp in self.low_throughput_rbps:
            Y[rbp][:] = np.nan
        # scan the data and update the interaction matrix
        for e in data.iterrows():
            rbp = e[1][0]
            rna = e[1][1]
            Y[rbp][rna] = 1
        # save the DataFrame to HDF5
        store = pd.io.pytables.HDFStore(self.output)
        store["matrix"] = Y
        store.close()

        if self.verbose:
            print("Done.\n")
            print("Interaction matrix saved to %s" % self.output)
            sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('interaction_list', metavar='interaction_list', type=str,
                        help="""File containing the interaction list. First column: proteins, second column: rnas (tab separated).""")
    parser.add_argument('output', metavar='output', type=str,
                        help="""File name of the HDF Store to save the intereaction matrix.""")
    parser.add_argument('--low', metavar='low_throughput_rbps', type=str, default=None,
                        help="""File containing a list (new line separated) of RBPs with low-throughput evidence only. Their missing interactions will be labeled with NaN.""")
    parser.add_argument('--quiet', dest='quiet', action='store_true', default=False,
                        help="""Do not print information at STDOUT.""")

    args = parser.parse_args()

    im = InteractionMatrix(interaction_list=args.interaction_list,
                        output=args.output,
                        verbose=(not args.quiet))
    im.prepare()
