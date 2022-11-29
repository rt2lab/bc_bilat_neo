#!/usr/bin/env python

import pandas as pd
import os

pd.options.display.max_columns = 200
rows, columns = os.popen('stty size', 'r').read().split()
pd.set_option('display.width', int(columns))
pd.set_option('display.max_columns', int(rows))

DECAY = 100

samtools_path = "/bioinfo/local/samtools/samtools"
bilat_path = "/Users/JudithAbecassis/Documents/PhD/bilat_ashp"
snapshot_directory = '{}/mutations_igv'.format(bilat_path)

server_script_name = 'get_mini_bam_mutations.sh'
igv_script_name = 'igv_script_mutations'
bam_table = pd.read_csv('superfreq_metadata_patient.tsv', sep='\t')

left_right_table = pd.read_csv('left_right_common_mutations.csv', sep='\t')
left_right_table_dedup = left_right_table.drop_duplicates(subset=['patient_id', 'mutation_id'])
left_right_table_dedup = left_right_table_dedup.assign(motif='left_right_intersection')
driver_table = pd.read_csv('driver_mutations.csv', sep='\t')
driver_table = driver_table.assign(motif='driver')
driver_table_dedup = driver_table.drop_duplicates(subset=['patient_id', 'mutation_id'])

mutation_table = pd.concat((left_right_table_dedup, driver_table_dedup), axis=0)

server_script = open(server_script_name, 'w')
igv_script = open(igv_script_name, 'w')

server_script.write('mkdir mini_bam_mutations\n')

for idx, row in mutation_table.iterrows():
    gene = row['inGene']
    pos = row['start']
    loc = 'chr{}:{}-{}'.format(row['chr'], pos - DECAY, pos + DECAY)
    relevant_bam = bam_table[bam_table.INDIVIDUAL.str.contains(str(row['patient_id']))]
    igv_script.write('new\nsnapshotDirectory {}\n'.format(snapshot_directory))
    for idx_bam, row_bam in relevant_bam.iterrows():
        server_script.write('{} view -b {} "{}" > mini_bam_mutations/{}_{}_chr{}_{}_{}.bam\n'.format(samtools_path, row_bam['BAM'], loc, row_bam['NAME'], gene, row['chr'], pos, row['motif']))
        server_script.write('{} index mini_bam_mutations/{}_{}_chr{}_{}_{}.bam\n\n'.format(samtools_path, row_bam['NAME'], gene, row['chr'], pos, row['motif']))

        igv_script.write('load {}/mini_bam_mutations/{}_{}_chr{}_{}_{}.bam\n'.format(bilat_path, row_bam['NAME'], gene, row['chr'], pos, row['motif']))
        igv_script.write('goto {}\n'.format(loc))
    igv_script.write('snapshot P{}_{}_chr{}_{}_{}.png\n\n'.format(row['patient_id'], gene, row['chr'], pos, row['motif']))


server_script.close()
igv_script.close()
