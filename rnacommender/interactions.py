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
