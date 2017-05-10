#!/usr/bin/env python3

from Bio import SeqIO
import csv
import tompytools

# read genbank file
genbank_file = 'output/gb/rbcl_58024.gb'
with open(genbank_file, 'r') as genbank_handle:
    gb_records = list(SeqIO.parse(genbank_file, 'gb'))

# rename to ID only
for record in gb_records:
    record.description = ''

# write output
fasta_file = 'output/fa/rbcl_58024.fa'
SeqIO.write(gb_records, fasta_file, 'fasta')

# write accession to taxonomy file
lines = [
    list(tompytools.flatten_list([x.id, x.annotations['taxonomy']]))
    for x in gb_records]

with open('output/taxonomy/genbank_results.txt', 'w') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerows(lines)
