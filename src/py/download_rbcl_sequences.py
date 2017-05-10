#!/usr/bin/env python3

import argparse
from Bio import Entrez
from Bio import SeqIO

# require email from cli
parser = argparse.ArgumentParser()
parser.add_argument('-e',
                    help='email address',
                    metavar='email',
                    type=str,
                    required=True)
args = parser.parse_args()

# identify myself to Entrez
Entrez.email = args.e

# spermatocytes: txid58024[Organism] 
search_term = ('rbcL[Title] '
               'AND txid58024[Organism]')

# initial search to get the number of hits and webenv
with Entrez.esearch(db='nucleotide',
                    term=search_term,
                    usehistory='y',
                    idtype='acc') as handle:
    search_results = Entrez.read(handle)
    number_of_hits = int(search_results['Count'])
    webenv = search_results['WebEnv']
    query_key = search_results['QueryKey']

print("%i hits" % number_of_hits)

# download and parse genes 1000 at a time
batch_size = 1000
with open('output/gb/rbcl_58024.gb', 'w') as genbank_handle:
    for start in range(0, number_of_hits, batch_size):
        end = min(number_of_hits, start+batch_size)
        print("Start:\t%s\nEnd:\t%s" % (start, end))

        fetch_handle = Entrez.efetch(
                db='nucleotide',
                idtype='acc',
                rettype='gb',
                retmode='text',
                webenv=webenv,
                query_key=query_key,
                retstart=start,
                retmax=batch_size)
        gb_records = SeqIO.parse(fetch_handle, 'gb')
        SeqIO.write(gb_records, genbank_handle, 'gb')
