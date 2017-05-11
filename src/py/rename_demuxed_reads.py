#!/usr/bin/env python3

from Bio import SeqIO
import tompytools
import gzip
import os
import re
import csv

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

# write sample mapping file
with open(os.path.join(outdir, "mapping.txt"), 'w') as f:
    f.write('#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n')
    lines = ['\t'.join([sample_names[x], expected_barcodes[x],
                        linker_primers[x], sample_names[x]])
             for x in sample_names]
    f.writelines('%s\n' % line for line in lines)

# parse sample names
basename_regex = r'^(?P<SN>2556-\d+)-.*'

# parse barcode
desc_regex = r'^.+\s.+\:(?P<BC>\w+\+\w+)$'

# find merged files
fastq_files = tompytools.find_all(
    ['_merged.fastq.gz'], 'output/trim_merge')

# collect modified fastq
modified_fastq = []

# open a handle for writing fasta output
with open(os.path.join(outdir, 'seqs.fna'), 'w') as outfile:
    # loop over fastq
    for fq in fastq_files:
        # get the sample name
        bn = os.path.basename(fq)
        nzgl_name = re.sub(basename_regex, r'\g<SN>', bn)
        orig_bc = expected_barcodes[nzgl_name]
        sample_name = sample_names[nzgl_name]
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
                new_bc = re.sub(desc_regex, r'\g<BC>', desc)
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

# write output to fasta
SeqIO.write(modified_fastq,
            'test/seqs.fna',
            'fasta')
