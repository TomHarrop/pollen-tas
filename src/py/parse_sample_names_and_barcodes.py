#!/usr/bin/env python3

import csv
from Bio.Seq import Seq

# parse barcode sequences
bc_seqs = {}
rc_bc_seqs = {}

with open('data/project_files/bc_name_and_seq.csv', 'r') as f:
    rows = csv.reader(f)
    next(rows)
    for row in rows:
        bc_seqs[row[0]] = Seq(row[1])

for key in bc_seqs:
    rc_bc_seqs[key] = bc_seqs[key].reverse_complement()

# parse sample name file
sample_names = {}
first_bc_rc = {}
second_bc = {}
with open('data/project_files/sample_id_and_bc.csv', 'r') as f:
    rows = csv.reader(f)
    next(rows)
    for row in rows:
        sample_names[row[0]] = row[1]
        first_bc_rc[row[0]] = row[2]
        second_bc[row[0]] = row[3]

# construct expected barcode
expected_barcode_seq = {}
for key in sample_names:
    expected_barcode_seq[key] = (
        "%s+%s" %
        (rc_bc_seqs[first_bc_rc[key]], bc_seqs[second_bc[key]]))

# write output
csvrows = [[x, sample_names[x], expected_barcode_seq[x]]
           for x in expected_barcode_seq]
with open('data/project_files/expected_barcodes.csv', 'w') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["nzgl_id", "sample_name", "expected_barcodes"])
    csvwriter.writerows(csvrows)


