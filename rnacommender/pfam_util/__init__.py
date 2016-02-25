from urllib2 import urlopen, URLError

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

def scan(seq):
    raise NotImplementedError("Implement wrapper for pfam_scan.pl")

def _dowload_seed_seqs(acc):

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

def download_seeds(dom_list,temp_folder,verbose):

