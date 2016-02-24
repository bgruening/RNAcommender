from __future__ import print_function

import argparse
import pandas as pd
import sys
import subprocess as sp
from urllib2 import urlopen, URLError
from uuid import uuid4
from os import mkdir,rmdir

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

"""
NOTES:
* move all the methods that interact with Pfam in a pfam util file (easier to update when pfam changes)
"""

class RBPVectorizer():
    def __init__(self,fasta_ref,fasta_sel,pfam_scan_ref=None,pfam_scan_sel=None):
        self.fasta_ref = fasta_ref
        self.fasta_sel = fasta_sel
        self.pfam_scan_ref = pfam_scan_ref
        self.pfam_scan_sel = pfam_scan_sel

    def pfam_scan(self):
        raise NotImplementedError("Implement wrapper for pfam_scan.pl")

    def _overlapping_domains(self):
        cols = ["seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E-value","significance","clan"]
        if self.pfam_scan_ref != self.pfam_scan_sel:
            data_ref = pd.read_table(self.pfam_scan_ref,
                sep="\s*",skip_blank_lines=True,skiprows=1,names=cols)
            doms_ref = set(a.split('.')[0] for a in data_ref["hmm_acc"])
            data_sel = pd.read_table(self.pfam_scan_sel,
                sep="\s*",skip_blank_lines=True,skiprows=1,names=cols)
            doms_sel = set(a.split('.')[0] for a in data_sel["hmm_acc"])
            return sorted(list(doms_ref & doms_sel))
        else:
            data = pd.read_table(self.pfam_scan_ref,
                sep="\s*",skip_blank_lines=True,skiprows=1,names=cols)
            data["hmm_acc"] = [a.split('.')[0] for a in data["hmm_acc"]]
            rbp_dom = data[["seq_id","hmm_acc"]].drop_duplicates()
            group = rbp_dom.groupby("hmm_acc").count()
            doms = (group[group>1]).dropna().index
            return sorted(list(doms))

    def _stockholm2fasta(self,stockholm):
        fasta = ""
        for line in stockholm.split("\n"):
            # comment line
            if line[0] == "#":
                continue
            # termination line
            elif line == "//":
                return fasta
            # alignment line
            else:
                name,align = line.split()
                seq = align.replace(".","")
                fasta += ">%s\n%s\n" % (name,seq)

    def _dowload_seed_seqs(self,acc):
        try:
            url = "http://pfam.xfam.org/family/%s/alignment/seed" % acc
            socket = urlopen(url)
        except URLError:
            raise URLError("Accession not recognized: %s (go to http://pfam.xfam.org/ for more details)." % acc)
        else:
            stockholm = socket.read()
            socket.close()
            fasta = self._stockholm2fasta(stockholm)
            return fasta

    def vectorize(self):
        # create a temporary hidden folder
        temp_fold = "temp_" + str(uuid4())
        mkdir(temp_fold)

        # if the file is not given as input run pfam scan
        if self.pfam_scan_ref is None or self.pfam_scan_sel is None:
            print("Scanning fasta files against the Pfam Database...", end=' ')
            sys.stdout.flush()
            self._pfam_scan() # raises NotImplementedError
            print("Done.\n")
            sys.stdout.flush()

        print("Determining domain list...", end=' ')
        sys.stdout.flush()
        # determine the accession numbers of the pfam domains
        dom_list = self._overlapping_domains()
        print("Done.\n")
        sys.stdout.flush()

        # download the alignment of the seeds from pfam and convert it to fasta
        seeds_fold = "%s/seeds" % temp_fold
        mkdir(seeds_fold)
        print("Downloading %i domain seeds from http://pfam.xfam.org/..." % len(dom_list), end=' ')
        sys.stdout.flush()
        for acc in dom_list:
            seed = self._dowload_seed_seqs(acc)
            if seed is not None:
                nf = open("%s/seeds/%s.fa" % (temp_fold,acc),"w")
                nf.write(seed)
                nf.close()
        print("Done.\n")
        sys.stdout.flush()

        # compile the models using SAM 3.5
        mod_fold = "%s/mod" % temp_fold
        mkdir(mod_fold)
        print("Building %i HMM models..." % len(dom_list))
        sys.stdout.flush()
        for i,acc in enumerate(dom_list):
            cmd = "buildmodel %s/%s -train %s/%s.fa -randseed 0" % (mod_fold,acc,seeds_fold,acc)
            sp.check_call(cmd,shell=True)
        print("Done.\n")
        sys.stdout.flush()


v = RBPVectorizer(fasta_ref="../examples/rbps.fa_HT",
    fasta_sel="../examples/rbps_HT.fa",
    pfam_scan_ref="../examples/pfam_scan_HT.txt",
    pfam_scan_sel="../examples/pfam_scan_HT.txt")
v.vectorize()

