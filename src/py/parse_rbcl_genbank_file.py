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

# wtf?
test_blank = [x for x in gb_records if x.id == '']

# write long
lines = []
for record in gb_records:
    for tax in record.annotations['taxonomy']:
        lines.append([record.annotations['organism'], record.id, tax])
with open('output/taxonomy/genbank_results.txt', 'w') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow(['species', 'accession', 'taxon'])
    csv_writer.writerows(lines)
