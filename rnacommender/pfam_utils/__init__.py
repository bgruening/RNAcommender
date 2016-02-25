import pandas as pd
from urllib2 import urlopen, URLError

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

def scan(seq):
    raise NotImplementedError("Implement wrapper for pfam_scan.pl")

def read_pfam_output(pfam_out_file):
    cols = ["seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E-value","significance","clan"]
    data = pd.read_table(pfam_out_file,
            sep="\s*",skip_blank_lines=True,skiprows=1,names=cols)
    return data

def download_seed_seqs(acc):

    def stockholm2fasta(stockholm):
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

    try:
        url = "http://pfam.xfam.org/family/%s/alignment/seed" % acc
        socket = urlopen(url)
    except URLError:
        raise URLError("Accession not recognized: %s (go to http://pfam.xfam.org/ for more details)." % acc)
    else:
        stockholm = socket.read()
        socket.close()
        fasta = stockholm2fasta(stockholm)
        return fasta
