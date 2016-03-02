import pandas as pd
from urllib import urlencode
from urllib2 import Request, urlopen, URLError
from httplib import HTTPException
import xml.etree.ElementTree as ET
from xml.etree.ElementTree import ParseError
from math import ceil
from time import sleep

import fasta_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"

def search_header():
    return "<seq id>        <alignment start>       <alignment end> <envelope start>        <envelope end>  <hmm acc>       <hmm name>      <type>  <hmm start>     <hmm end>       <hmm length>    <bit score>     <E-value>       <significance>  <clan>\n"

def sequence_search(seq_id,seq):

    def add_spaces(text,mul=8):
        l = len(text)
        next_mul = int(ceil(l / mul)+1) * mul
        offset = next_mul - l
        if offset == 0:
            offset = 8
        return text + " "*offset

    url = "http://pfam.xfam.org/search/sequence"
    params = {'seq' : seq,
            'evalue' : '1.0',
            'output' : 'xml' }
    data = urlencode(params)
    req = Request(url, data)
    # give some time to Pfam to elaborate the request
    sleep(0.5)
    while urlopen(req).getcode() != 200:
        sleep(0.5)
    response = urlopen(req)
    xml = response.read()
    root = ET.fromstring(xml)
    result_url = root[0][1].text
    # wait for Pfam to compute the results
    sleep(4)
    while urlopen(result_url).getcode() != 200:
        sleep(1)
    socket = urlopen(result_url)
    result_xml = socket.read()
    root = ET.fromstring(result_xml)
    matches = root[0][0][0][0][:]
    ret = ""
    for match in matches:
        for location in match:
            ret += add_spaces(seq_id)
            ret += add_spaces(location.attrib['ali_start'])
            ret += add_spaces(location.attrib['ali_end'])
            ret += add_spaces(location.attrib['start'])
            ret += add_spaces(location.attrib['end'])
            ret += add_spaces(match.attrib['accession'])
            ret += add_spaces(match.attrib['id'])
            ret += add_spaces(match.attrib['class'])
            ret += add_spaces(location.attrib['hmm_start'])
            ret += add_spaces(location.attrib['hmm_end'])
            ret += add_spaces("None")
            ret += add_spaces(location.attrib['bitscore'])
            ret += add_spaces(location.attrib['evalue'])
            ret += add_spaces(location.attrib['significant'])
            ret += "None\n"
    return ret

def read_pfam_output(pfam_out_file):
    cols = ["seq_id","alignment_start","alignment_end","envelope_start","envelope_end","hmm_acc","hmm_name","type","hmm_start","hmm_end","hmm_length","bit_score","E-value","significance","clan"]
    data = pd.read_table(pfam_out_file,
            sep="\s*",skip_blank_lines=True,skiprows=1,names=cols,engine='python')
    return data

def download_seed_seqs(acc):
    try:
        url = "http://pfam.xfam.org/family/%s/alignment/seed" % acc
        socket = urlopen(url)
    except URLError:
        raise URLError("Accession not recognized: %s (go to http://pfam.xfam.org/ for more details)." % acc)
    else:
        stockholm = socket.read()
        socket.close()
        fasta = fasta_utils.stockholm2fasta(stockholm)
        return fasta
