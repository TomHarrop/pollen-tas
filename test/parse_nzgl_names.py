#!/usr/bin/env python3

from Bio import SeqIO
import re
import gzip
import os
import csv

# read expected_barcodes into dicts
sample_name = {}
expected_barcodes = {}
with open('data/expected_barcodes.csv', 'r') as f:
    rows = csv.reader(f)
    for row in rows:
        print(row)
        sample_name[row[0]] = row[1]
        expected_barcodes[row[0]] = row[2]

# parse the nzgl id from the filename
filename = 'data/NZGL02556/2556-25-20-01_S25_L001_R1_001.fastq.gz'
filename_grep = r'^(?P<NN>2556-\d+)-.*'
nzgl_name = re.sub(filename_grep,
       r'\g<NN>',
       os.path.basename(filename))

# set up expected values
exp_bc = expected_barcodes[nzgl_name]
sn = sample_name[nzgl_name]

# barcode grep
bc_grep = r'^.+\s.+\:(?P<BC>\w+\+\w+)$'

# read the fastq
with gzip.open(filename, 'rt') as handle:
    seq = SeqIO.parse(handle, 'fastq-sanger')
    records = [x for x in seq]
    # modify the records
    i=0
    for record in records:
        i+=1
        # actual barcode
        read_bc = re.sub(bc_grep, r'\g<BC>', record.description)
        # number of differences between read_bc and exp_bc
        bc_diff = sum(c1!=c2
                      for c1,c2
                      in zip(exp_bc, read_bc))
        # construct a new id
        new_id = (
            "%s_%i exp_bc=%s read_bc=%s bc_diff=%i"
            % (sn, i, exp_bc, read_bc, bc_diff))
        # to do: replace ID, add to record, write to fasta

