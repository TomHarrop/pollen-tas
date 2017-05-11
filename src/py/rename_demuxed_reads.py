#!/usr/bin/env python3

from Bio import SeqIO
import tompytools
import gzip
import os
import re
import csv

# parse expected_barcodes and nzgl names from project files
sample_names = {}
expected_barcodes = {}
with open('data/project_files/expected_barcodes.csv', 'r') as f:
    rows = csv.reader(f)
    for row in rows:
        sample_names[row[0]] = row[1]
        expected_barcodes[row[0]] = row[2]

# parse sample names
basename_regex = r'(?P<SN>^.+_S\d+)_.*'

# parse barcode
desc_regex = r'.*\:(?P<BC>.+)$'

# find merged files
fastq_files = tompytools.find_all(
    ['_merged.fastq.gz'], 'output/trim_merge')

# collect modified fastq
modified_fastq = []

# collect sample name and barcodes for mapping file
barcode_to_sample = {}
sample_to_barcode = {}

# TODO:
# read the actual barcode and desired sample name from andrew's xlsx. calculate
# the "distance" between the barcode and the theoretical barcode in the loop

# loop over fastq
for fq in fastq_files:
    # get the sample name
    bn = os.path.basename(fq)
    sample_name = re.sub(basename_regex, r'\g<SN>', bn)
    # unzip and read fq
    with gzip.open(fq, 'rt') as handle:
        fastq = SeqIO.parse(handle, 'fastq-sanger')
        records = [x for x in fastq]
        # modify records
        i=0
        for record in records:
            i+=1
            # get the barcode
            desc = record.description
            barcode = re.sub(desc_regex, r'\g<BC>', desc)
            # collect barcode for mapping
            barcode_to_sample[barcode] = sample_name
            sample_to_barcode[sample_name] = barcode
            # write the new id string
            new_id_string = (
                "%s_%i %s orig_bc=%s new_bc=%s bc_diffs=0" %
                (sample_name, i, record.id, barcode, barcode))
            record.id = new_id_string
            record.name = ''
            record.description = ''
            # collect result
            modified_fastq.append(record)

# write output to fasta
SeqIO.write(modified_fastq,
            'test/seqs.fna',
            'fasta')
