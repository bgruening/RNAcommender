
__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

def import_fasta(fasta_file):
    """
    Import a fasta file as a dictionary, k:v, where k is the
    sequence name, and v is the sequence.
    """
    dic = {}
    f = open(fasta_file)
    fasta = f.read().strip()
    f.close()
    for a in fasta.split('>'):
        k = a.split('\n')[0]
        v = ''.join(a.split('\n')[1:])
        if k != '':
            dic[k] = v
    return dic

def stockholm2fasta(stockholm):
    """
    Convert alignment in stockholm format to
    fasta format.
    """
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