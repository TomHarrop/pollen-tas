#!/usr/bin/env python3

from Bio import SeqIO
import tompytools
import gzip
import os
import re
import csv


# barcode sanitising wrapper
def sanitise_barcode(x):
    return(re.sub('\+', '', x))


# name sanitising wrapper
def sanitise_name(x):
    return(re.sub('[^\w]|_', '.', x))


# output:
outdir = "output/qiime/split_libraries"

# parse expected_barcodes and nzgl names from project files
sample_names = {}
expected_barcodes = {}
with open('data/project_files/expected_barcodes.csv', 'r') as f:
    rows = csv.reader(f)
    next(rows)
    for row in rows:
        sample_names[row[0]] = row[1]
        expected_barcodes[row[0]] = row[2]

# read primers
linker_primers = {}
with open('data/project_files/linker_primers.csv', 'r') as f:
    rows = csv.reader(f)
    for row in rows:
        linker_primers[row[0]] = row[1]

# remove plus symbol from barcode
sanitised_barcodes = {}
for x in expected_barcodes:
    sanitised_barcodes[x] = sanitise_barcode(expected_barcodes[x])

# remove disallowed characters from names
sanitised_names = {}
for x in sample_names:
    sanitised_names[x] = sanitise_name(sample_names[x])

# write sample mapping file
with open(os.path.join(outdir, "mapping.txt"), 'w') as f:
    f.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n')
    lines = ['\t'.join([sanitised_names[x], sanitised_barcodes[x],
                        linker_primers[x], sample_names[x]])
             for x in sample_names]
    f.writelines('%s\n' % line for line in lines)

# parse sample names
basename_regex = r'^D2JR8-(?P<SN>2556-\d+)-.*'

# parse barcode
desc_regex = r'^.+\s.+\:(?P<BC>\w+\+\w+)$'

# find R1 files
fastq_files = tompytools.find_all(
    ['R1_001.fastq.gz'], 'data/NZGL02556')

# open a handle for writing fasta output
with open(os.path.join(outdir, 'seqs.fna'), 'w') as outfile:
    # loop over fastq
    for fq in fastq_files:
        # get the sample name
        bn = os.path.basename(fq)
        nzgl_name = re.sub(basename_regex, r'\g<SN>', bn)
        orig_bc = sanitised_barcodes[nzgl_name]
        sample_name = sanitised_names[nzgl_name]
        # unzip and read fq
        with gzip.open(fq, 'rt') as handle:
            fastq = SeqIO.parse(handle, 'fastq-sanger')
            records = [x for x in fastq]
            # modify records
            i = 0
            for record in records:
                i += 1
                # get the barcode
                desc = record.description
                new_bc = sanitise_barcode(re.sub(desc_regex, r'\g<BC>', desc))
                # number of differences between read_bc and exp_bc
                bc_diff = sum(c1 != c2
                              for c1, c2
                              in zip(orig_bc, new_bc))
                # write the new id string
                new_id_string = (
                    "%s_%i %s orig_bc=%s new_bc=%s bc_diffs=%i" %
                    (sample_name, i, record.id, orig_bc, new_bc, bc_diff))
                record.id = new_id_string
                record.name = ''
                record.description = ''
                # write result
                if bc_diff == 0:
                    SeqIO.write(record, outfile, 'fasta')
