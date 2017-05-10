#!/usr/bin/env python3

from Bio import SeqIO

records = [x for x in SeqIO.parse('previous_data/seq45R1.fasta',
                                  'fasta')]
for record in records:
    new_id_string = ("Test %s orig_bc=TGACCATTTCAC "
                     "new_bc=TGACCATTTCAC bc_diffs=0" % record.id)
    record.id = new_id_string
    record.name = ''
    record.description = ''

SeqIO.write(records,
            'output/qiime/split_libraries/seqs.fna',
            'fasta')
