from __future__ import print_function

import argparse
import pandas as pd
import sys
import subprocess as sp
from uuid import uuid4
from os import mkdir,rmdir

import pfam_interface

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class RBPVectorizer():
    def __init__(self,fasta_ref,fasta_sel,pfam_scan_ref=None,pfam_scan_sel=None):
        self.fasta_ref = fasta_ref
        self.fasta_sel = fasta_sel
        self.pfam_scan_ref = pfam_scan_ref
        self.pfam_scan_sel = pfam_scan_sel

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

    def _prepare_domains(self,dom_list,dom_ref_fold,dom_sel_fold):

        def import_fasta(fasta_file):
            dic = {}
            f = open(fasta_file)
            fasta = f.read()
            f.close()
            for a in fasta.split('>'):
                k = a.split('\n')[0]
                v = ''.join(a.split('\n')[1:])
                dic[k] = v
            return dic

        def prepare_domains(fasta_dic,dom_list,pfam_scan,out_folder):
            out_file_dic = {}
            for acc in dom_list:
                out_file_dic[acc] = open("%s/%s.fa" % (out_folder,acc), "w")

            f = open(pfam_scan)
            f.readline()
            for line in f:
                split = line.split()
                rbp = split[0]
                start = int(split[3])
                stop = int(split[4])
                acc = split[5].split('.')[0]
                if acc in out_file_dic.keys():
                    out_file_dic[acc].write(">%s:%i-%i\n%s\n" % (rbp,start,stop,fasta_dic[rbp][start:stop]))
            f.close()

            for acc in dom_list:
                out_file_dic[acc].close()

        fasta_ref = import_fasta(self.fasta_ref)
        prepare_domains(fasta_ref,dom_list,self.pfam_scan_ref,dom_ref_fold)

        fasta_sel = import_fasta(self.fasta_sel)
        prepare_domains(fasta_sel,dom_list,self.pfam_scan_sel,dom_sel_fold)

    def vectorize(self):
        # create a temporary hidden folder
        temp_fold = "temp_" + str(uuid4())
        mkdir(temp_fold)

        # run pfam scan (if needed)
        if self.pfam_scan_ref is None or self.pfam_scan_sel is None:
            raise NotImplementedError()

        print("Determining domain list...", end=' ')
        sys.stdout.flush()
        # determine the accession numbers of the pfam domains
        dom_list = self._overlapping_domains()
        print("Done.\n")
        sys.stdout.flush()

        #prepare fasta file with the sequence of the domains
        dom_ref_fold = "%s/domains_ref" % temp_fold
        mkdir(dom_ref_fold)
        dom_sel_fold = "%s/domains_sel" % temp_fold
        mkdir(dom_sel_fold)
        print("Preparing fasta files for %i domains..." % len(dom_list), end=' ')
        sys.stdout.flush()
        self._prepare_domains(dom_list,dom_ref_fold,dom_sel_fold)
        print("Done.\n")
        sys.stdout.flush()

        # download the alignment of the seeds from pfam and convert it to fasta
        seeds_fold = "%s/seeds" % temp_fold
        mkdir(seeds_fold)
        print("Downloading %i domain seeds from http://pfam.xfam.org/..." % len(dom_list), end=' ')
        sys.stdout.flush()
        for acc in dom_list:
            seed = pfam_interface.dowload_seed_seqs(acc)
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


v = RBPVectorizer(fasta_ref="../examples/rbps_HT.fa",
    fasta_sel="../examples/rbps_HT.fa",
    pfam_scan_ref="../examples/pfam_scan_HT.txt",
    pfam_scan_sel="../examples/pfam_scan_HT.txt")
v.vectorize()

