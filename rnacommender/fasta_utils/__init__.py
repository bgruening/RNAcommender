def import_fasta(fasta_file):
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