from __future__ import print_function

import argparse
import numpy as np
import pandas as pd
import sys
import subprocess as sp
import uuid
import re
from os import mkdir,rmdir

import pfam_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


class RBPVectorizer():
    def __init__(self,fasta_ref,fasta_sel,output,pfam_scan_ref=None,pfam_scan_sel=None,include_all_sel=False,verbose=True):
        self.fasta_ref = fasta_ref
        self.fasta_sel = fasta_sel
        self.output = output
        self.pfam_scan_ref = pfam_scan_ref
        self.pfam_scan_sel = pfam_scan_sel
        self.include_all_sel = include_all_sel
        self.verbose = verbose

        # self.temp_fold = "temp_" + str(uuid.uuid4())
        self._temp_fold = "temp_38267b24-0c13-47cb-b66c-5e71184d52d9"
        self._dom_ref_fold = "%s/domains_ref" % self._temp_fold
        self._dom_sel_fold = "%s/domains_sel" % self._temp_fold
        self._seeds_fold = "%s/seeds" % self._temp_fold
        self._mod_fold = "%s/mod" % self._temp_fold
        self._fisher_ref_fold = "%s/fisher_scores_ref" % self._temp_fold
        self._fisher_sel_fold = "%s/fisher_scores_sel" % self._temp_fold

    def _overlapping_domains(self):
        if self.verbose:
            print("Determining domain list...", end=' ')
            sys.stdout.flush()

        if self.pfam_scan_ref != self.pfam_scan_sel:
            data_ref = pfam_utils.read_pfam_output(self.pfam_scan_ref)
            dom_ref = set(a.split('.')[0] for a in data_ref["hmm_acc"])
            data_sel = pfam_utils.read_pfam_output(self.pfam_scan_sel)
            dom_sel = set(a.split('.')[0] for a in data_sel["hmm_acc"])
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

        if self.verbose:
            print("Preparing fasta files with domain sequences...", end=' ')
            sys.stdout.flush()

        mkdir(self._dom_ref_fold)
        fasta_ref = import_fasta(self.fasta_ref)
        prepare_domains(fasta_ref,dom_list,self.pfam_scan_ref,self._dom_ref_fold)

        mkdir(self._dom_sel_fold)
        fasta_sel = import_fasta(self.fasta_sel)
        prepare_domains(fasta_sel,dom_list,self.pfam_scan_sel,self._dom_sel_fold)

        if self.verbose:
            print("Done.\n")
            sys.stdout.flush()

    def _download_seeds(self,dom_list):
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

        def process_seg(e):
            seg = e.split()
            c = seg[0].split(':')[0]
            m = map(float,seg[3:])
            return c,m

        def read_sam_file(samfile):
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
        # create a temporary hidden folder
        mkdir(self._temp_fold)

        # pfam scan is not implemented yet
        if self.pfam_scan_ref is None or self.pfam_scan_sel is None:
            raise NotImplementedError("Run pfam_scan from: http://pfam.xfam.org/search#tabview=tab1")

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
        compute the empirical kernel map
        return self._ekm(dom_list)



v = RBPVectorizer(fasta_ref="../examples/rbps_HT.fa",
    fasta_sel="../examples/rbps_HT.fa",
    pfam_scan_ref="../examples/pfam_scan_HT.txt",
    pfam_scan_sel="../examples/pfam_scan_HT.txt",
    include_all_sel=True,
    output="prova.h5")
a = v.vectorize()

