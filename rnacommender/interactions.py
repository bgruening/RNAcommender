#!/usr/bin/env python

from __future__ import print_function
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
    def __init__(self,interaction_list,low_throughput_rbps=None,verbose=True):
        """
        Params
        ------
        interaction_list : str
            File containing the interaction list. Comments and header
            must start with #.

        low_throughput_rbps : str (default : None)
            File containing a list (new line separated) of RBPs
            with low-throughput evidence only. Their missing interactions
            will be labeled with NaN.

        verbose : bool (default : True)
            Print information to STDOUT.
        """
        self.interaction_list = interaction_list
        if low_throughput_rbps is not None:
            f = open(low_throughput_rbps)
            self.low_throughput_rbps = f.read().strip().split()
            f.close()
        else:
            self.low_throughput_rbps = []
        self.verbose = verbose

    def prepare(self):
        #read interaction list
        data = pd.read_table(self.interaction_list,sep='\t',comment='#',header=None)
        # determine the names of the RBPs and the RNAs
        rbps = data[0].unique()
        rnas = data[1].unique()
        # Define the DataFrame with all the intarction at 0
        n = rbps.shape[0]
        m = rnas.shape[0]
        df = pd.DataFrame(np.zeros((m,n)),columns=rbps,index=rnas)
        # Set NaN the interactions for the low-throughput RBPs
        for rbp in self.low_throughput_rbps:
            df[rbp][:] = np.nan

        # check df.iterrows()

i = InteractionMatrix(interaction_list="../examples/interactions_small.txt")
a = i.prepare()
