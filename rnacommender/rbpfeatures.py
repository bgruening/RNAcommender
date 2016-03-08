#!/usr/bin/env python

from __future__ import print_function
import argparse
import numpy as np
import pandas as pd
import sys
import subprocess as sp
import uuid
import re
from os import mkdir
from shutil import rmtree

import pfam_utils
import fasta_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class RBPVectorizer():
    """Computes the RBP features"""
    def __init__(self,fasta_ref,fasta_sel,output,include_all_sel=False,verbose=True):
        """
        Parameters
        ----------
        fasta_ref : str
            Fasta file containing the reference sequences. The similarity will computed
            against the reference sequences.

        fasta_sel : str
            Fasta file containing the selected sequences. The similarity will computed
            for the selected sequences. (This might be the same file as fasta_ref).

        output : str
            Name of the output file. The output file is an HDF Store containin a
            pandas DataFrame, where the columns are the selected sequence names and the rows
            are the reference sequence names

        include_all_sel : bool (default: False)
            Includes all the selected sequences even when they have zero similarity with all the
            reference sequences. If a sequence is both in the reference and selected set, and it
            has zero similarity with all the reference proteins except itself it will included only if
            this flag is set.

        verbose : bool (default : True)
            Print information to STDOUT
        """
        self.fasta_ref = fasta_ref
        self.fasta_sel = fasta_sel
        self.output = output
        self.pfam_scan_ref = None
        self.pfam_scan_sel = None
        self.include_all_sel = include_all_sel
        self.verbose = verbose

        self._temp_fold = "temp_" + str(uuid.uuid4())
        self._dom_ref_fold = "%s/domains_ref" % self._temp_fold
        self._dom_sel_fold = "%s/domains_sel" % self._temp_fold
        self._seeds_fold = "%s/seeds" % self._temp_fold
        self._mod_fold = "%s/mod" % self._temp_fold
        self._fisher_ref_fold = "%s/fisher_scores_ref" % self._temp_fold
        self._fisher_sel_fold = "%s/fisher_scores_sel" % self._temp_fold

    def _pfam_scan(self):
        """Scan the sequences in fasta_ref and fasta_sel against the Pfam database"""
        if self.verbose:
            print("Scanning RBP sequences against Pfam...")
            sys.stdout.flush()

        if self.fasta_ref != self.fasta_sel:
            self.pfam_scan_ref = "%s/pfam_scan_ref.txt" % self._temp_fold
            self.pfam_scan_sel = "%s/pfam_scan_sel.txt" % self._temp_fold
            nf_ref = open(self.pfam_scan_ref, "w")
            nf_ref.write(pfam_utils.search_header())
            nf_sel = open(self.pfam_scan_sel, "w")
            nf_sel.write(pfam_utils.search_header())

            fasta_ref = fasta_utils.import_fasta(self.fasta_ref)
            fasta_sel = fasta_utils.import_fasta(self.fasta_sel)

            to_scan = sorted(set(fasta_ref.keys() + fasta_sel.keys()))
            for rbp in to_scan:
                if self.verbose:
                    print(rbp, end=', ')
                    sys.stdout.flush()

                if rbp in fasta_ref.keys() and rbp in fasta_sel.keys():
                    try:
                        assert fasta_ref[rbp] == fasta_sel[rbp]
                    except AssertionError:
                        print('%s: sequence mismatch between ref and sel')
                    seq = fasta_ref[rbp]
                    text = pfam_utils.sequence_search(rbp,seq)
                    nf_ref.write(text)
                    nf_sel.write(text)
                elif rbp in fasta_ref.keys():
                    seq = fasta_ref[rbp]
                    text = pfam_utils.sequence_search(rbp,seq)
                    nf_ref.write(text)
                else:
                    seq = fasta_sel[rbp]
                    text = pfam_utils.sequence_search(rbp,seq)
                    nf_sel.write(text)

            nf_ref.close()
            nf_sel.close()
        else:
            pfam_scan_file = "%s/pfam_scan.txt" % self._temp_fold
            self.pfam_scan_ref = pfam_scan_file
            self.pfam_scan_sel = pfam_scan_file

            nf = open(pfam_scan_file, "w")
            nf.write(pfam_utils.search_header())

            fasta = fasta_utils.import_fasta(self.fasta_ref)

            for rbp in sorted(fasta.keys()):
                if self.verbose:
                    print(rbp, end=', ')
                    sys.stdout.flush()
                seq = fasta[rbp]
                text = pfam_utils.sequence_search(rbp,seq)
                nf.write(text)

            nf.close()


        if self.verbose:
                print("Done.\n")
                sys.stdout.flush()

    def _overlapping_domains(self):
        """Compute the set of overlapping domains between the proteins in fasta_ref and fasta_sel"""
        if self.verbose:
            print("Determining domain list...", end=' ')
            sys.stdout.flush()

        if self.pfam_scan_ref != self.pfam_scan_sel:
            data_ref = pfam_utils.read_pfam_output(self.pfam_scan_ref)
            doms_ref = set(a.split('.')[0] for a in data_ref["hmm_acc"])
            data_sel = pfam_utils.read_pfam_output(self.pfam_scan_sel)
            doms_sel = set(a.split('.')[0] for a in data_sel["hmm_acc"])
            dom_list = sorted(list(doms_ref & doms_sel))
        else:
            data = pfam_utils.read_pfam_output(self.pfam_scan_ref)
            data["hmm_acc"] = [a.split('.')[0] for a in data["hmm_acc"]]
            rbp_dom = data[["seq_id","hmm_acc"]].drop_duplicates()
            group = rbp_dom.groupby("hmm_acc").count()
            doms = (group[group>1]).dropna().index
            dom_list = sorted(list(doms))

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

        return dom_list

    def _prepare_domains(self,dom_list):
        """Select domain sequences from the entire protein sequences"""

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

        if self.verbose:
            print("Preparing fasta files with domain sequences...", end=' ')
            sys.stdout.flush()

        mkdir(self._dom_ref_fold)
        fasta_ref = fasta_utils.import_fasta(self.fasta_ref)
        prepare_domains(fasta_ref,dom_list,self.pfam_scan_ref,self._dom_ref_fold)

        mkdir(self._dom_sel_fold)
        fasta_sel = fasta_utils.import_fasta(self.fasta_sel)
        prepare_domains(fasta_sel,dom_list,self.pfam_scan_sel,self._dom_sel_fold)

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

    def _download_seeds(self,dom_list):
        """Download seed sequences for the needed domains"""
        if self.verbose:
            print("Downloading domain seeds from http://pfam.xfam.org/...", end=' ')
            sys.stdout.flush()

        mkdir(self._seeds_fold)

        for acc in dom_list:
            seed = pfam_utils.download_seed_seqs(acc)
            if seed is not None:
                nf = open("%s/%s.fa" % (self._seeds_fold,acc),"w")
                nf.write(seed)
                nf.close()

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

    def _build_models(self,dom_list):
        """Wrapper for SAM 3.5 buildmodel"""
        if self.verbose:
            print("Building HMM models...")
            sys.stdout.flush()

        mkdir(self._mod_fold)

        for acc in dom_list:
            cmd = "buildmodel %s/%s -train %s/%s.fa -randseed 0" % (self._mod_fold,acc,self._seeds_fold,acc)
            sp.check_call(cmd,shell=True)

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

    def _compute_fisher_scores(self,dom_list):
        """Wrapper for SAM 3.5 get_fisher_scores"""
        def get_fisher_scores(dom_list,mod_fold,dom_fold,fisher_fold):
            for acc in dom_list:
                cmd = "get_fisher_scores run -i %s/%s.mod -db %s/%s.fa" % (mod_fold,acc,dom_fold,acc)
                fisher = sp.check_output(cmd,shell=True)
                nf = open("%s/%s.txt" % (fisher_fold,acc),"w")
                nf.write(fisher)
                nf.close()

        if self.verbose:
            print("Computing Fisher scores...")
            sys.stdout.flush()

        mkdir(self._fisher_ref_fold)
        get_fisher_scores(dom_list,self._mod_fold,self._dom_ref_fold,self._fisher_ref_fold)

        mkdir(self._fisher_sel_fold)
        get_fisher_scores(dom_list,self._mod_fold,self._dom_sel_fold,self._fisher_sel_fold)

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

    def _ekm(self,dom_list):
        """Compute the empirical kernel map from the Fisher scores"""
        def process_seg(e):
            """Process one segment of a SAM 3.5 get_fisher_scores output file"""
            seg = e.split()
            c = seg[0].split(':')[0]
            m = map(float,seg[3:])
            return c,m

        def read_sam_file(samfile):
            """Read a SAM 3.5 get_fisher_scores output file"""
            f = open(samfile)
            data = f.read()
            f.close()

            columns = []
            M = []
            split = re.split(">A ",data)[1:]
            for e in split:
                c,m = process_seg(e)
                columns.append(c)
                M.append(m)

            M = np.matrix(M)
            df = pd.DataFrame(data=M.T, columns=columns)
            return df

        def dom_features(fisher_fold,dom_list,names=None):
            """Compute the features with respect to a domain type"""
            dfs = []
            for acc in dom_list:
                df = read_sam_file("%s/%s.txt" % (fisher_fold,acc))
                df = df.groupby(df.columns,axis=1).mean()
                dfs.append(df)

            con = pd.concat(dfs,ignore_index=True)

            if names is not None:
                add = sorted(list(set(names)-set(con.columns)))
                fil = sorted(list(set(names)-set(add)))
                con = con[fil]
                for c in add:
                    con[c] = np.zeros(len(con.index),dtype='float64')
                con = con[names]

            con = con.fillna(0.0)
            return con

        def seq_names(fasta_file):
            """Get sequence names from fasta file"""
            names = []
            f = open(fasta_file)
            fasta = f.read()
            f.close()
            for a in fasta.split('>'):
                names.append(a.split('\n')[0])
            return [a for a in names if a != '']

        if self.verbose:
            print("Computing Fisher scores...", end=' ')
            sys.stdout.flush()

        ref_names = seq_names(self.fasta_ref)
        ref = dom_features(self._fisher_ref_fold,dom_list,names=ref_names)
        ekm_ref = ref.T.dot(ref)
        ekm_ref.index = ekm_ref.columns

        if self.include_all_sel:
            sel_names = seq_names(self.fasta_sel)
            sel = dom_features(self._fisher_sel_fold,dom_list, names=sel_names)
        else:
            sel = dom_features(self._fisher_sel_fold,dom_list)
        ekm_sel = sel.T.dot(sel)
        ekm_sel.index = ekm_sel.columns

        ekm = ref.T.dot(sel)

        for rs in ekm.columns:
            for rr in ekm.index:
                if ekm_ref[rr][rr] > 0 and ekm_sel[rs][rs] > 0:
                    ekm[rs][rr] /= np.sqrt(ekm_ref[rr][rr]*ekm_sel[rs][rs])

        if self.include_all_sel:
            # needed if a protein that is in ref and sel has no domains from dom_list
            for rs in ekm.columns:
                for rr in ekm.index:
                    if rr == rs:
                        ekm[rs][rr] = 1.0

        store = pd.io.pytables.HDFStore(self.output)
        store["features"] = ekm
        store.close()

        if self.verbose:
            print("Done.\n")
            print("RBP features saved in %s" % self.output)
            sys.stdout.flush()

    def vectorize(self):
        """Produce the RBP features"""
        # create a temporary hidden folder
        mkdir(self._temp_fold)
        # scan the RBP sequences against Pfam
        self._pfam_scan()
        # determine the accession numbers of the pfam domains needed for computing the features
        dom_list = self._overlapping_domains()
        #prepare fasta file with the sequence of the domains
        self._prepare_domains(dom_list)
        # download the alignment of the seeds from pfam and convert it to fasta
        self._download_seeds(dom_list)
        # compile the models using SAM 3.5
        self._build_models(dom_list)
        # compute fisher scores using SAM 3.5
        self._compute_fisher_scores(dom_list)
        #compute the empirical kernel map
        self._ekm(dom_list)
        # create a temporary hidden folder
        rmtree(self._temp_fold)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('fasta_ref', metavar='fasta_ref', type=str,
                        help="""Fasta file containing the reference RBP sequences.""")
    parser.add_argument('fasta_sel', metavar='fasta_sel', type=str,
                        help="""Fasta file containing the selected RBP sequences.""")
    parser.add_argument('output', metavar='output', type=str,
                        help="""File name of the HDF Store containing the RBP features.""")
    parser.add_argument('--all-sel', dest='all_sel', action='store_true', default=False,
                        help="""Return one vector for each selected RBP (even if the similarity is null with all the reference RBPs).""")
    parser.add_argument('--quiet', dest='quiet', action='store_true', default=False,
                        help="""Do not print information at STDOUT""")

    args = parser.parse_args()

    v = RBPVectorizer(fasta_ref=args.fasta_ref,
                    fasta_sel=args.fasta_sel,
                    include_all_sel=args.all_sel,
                    output=args.output,
                    verbose=(not args.quiet))
    v.vectorize()
