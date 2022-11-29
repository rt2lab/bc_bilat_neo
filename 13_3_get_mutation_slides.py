#!/usr/bin/env python

import pandas as pd
import os

pd.options.display.max_columns = 200
rows, columns = os.popen('stty size', 'r').read().split()
pd.set_option('display.width', int(columns))
pd.set_option('display.max_columns', int(rows))

DECAY = 100

bilat_path = "./PhD/bilat_ashp"
slides_filename = "igv_left_right.tex"

left_right_table = pd.read_csv('left_right_common_mutations_annot.csv', sep='\t')
#left_right_table_dedup = left_right_table.drop_duplicates(subset=['patient_id', 'mutation_id'])
left_right_table_dedup = left_right_table[~pd.isna(left_right_table.annotation_judith)]
left_right_table_dedup = left_right_table_dedup.assign(motif='left_right_intersection')

useful_cols = ['new_sample', 'chr', 'start', 'end', 'reference', 'variant', 'inGene', 'effect', 'cov', 'ref', 'var', 'dbMAF', 'ExAC_AF']

with open('{}/{}'.format(bilat_path, slides_filename), 'w') as slides:
    slides.write('\\documentclass[10pt]{beamer}\n')
    slides.write('\n')
    slides.write('\\usetheme[progressbar=frametitle]{metropolis}\n')
    slides.write('\\usepackage{appendixnumberbeamer}\n')
    slides.write('\\usepackage{lipsum}\n')
    slides.write('\\usepackage{amsmath}\n')
    slides.write('\\usepackage{amssymb}\n')
    slides.write('\\usepackage{booktabs}\n')
    slides.write('\\usepackage{hyperref}\n')
    slides.write('\\usepackage[scale=2]{ccicons}\n')
    slides.write('\\usepackage{pgfplots}\n')
    slides.write('\\usepgfplotslibrary{dateplot}\n')
    slides.write('\\usepackage{dirtytalk}\n')
    slides.write('\\usepackage{xspace}\n')
    slides.write('\\newcommand{\\themename}{\\textbf{\\textsc{metropolis}}\\xspace}\n')
    slides.write('\\definecolor{mpigreen}{HTML}{007977}\n')
    slides.write('\\setbeamercolor{frametitle}{bg=mpigreen}\n')
    slides.write('\n')
    slides.write('\\title{Inspection of left-right mutations}\n')
    slides.write('\n')
    slides.write('\\author{Judith}\n')
    slides.write('\\institute{RT2 Lab - Institut Curie}\n')
    slides.write('\n')
    slides.write('\\begin{document}\n')
    slides.write('\n')
    slides.write('\n')
    slides.write('\\maketitle\n')
    slides.write('\n')

    for idx, row in left_right_table_dedup.iterrows():
        gene = row['inGene']
        pos = row['start']
        loc = 'chr{}:{}-{}'.format(row['chr'], pos - DECAY, pos + DECAY)
        annot = row['annotation_judith'].replace('%', '\%')

        igv_filename = 'P{}_{}_chr{}_{}_{}.png'.format(row['patient_id'], gene, row['chr'], pos, row['motif'])

        slides.write('\\begin{frame}\n')
        slides.write('\\begin{figure}\n')
        slides.write('\\centering\n')
        slides.write('\\includegraphics[height=0.7\\textheight]{{mutations_igv/{}}}\n'.format(igv_filename))
        slides.write('\\end{figure}\n')

        table = left_right_table[(left_right_table.chr==row['chr'])&(left_right_table.start==pos)][useful_cols].to_latex(index=False)
        slides.write('\\Tiny{{{}}}\n'.format(table))

        slides.write('\n\n')
        slides.write('\\normalsize{{\\textcolor{mpigreen}{' + annot + '}}}\n')
        slides.write('\\end{frame}\n')

    slides.write('\n')
    slides.write('\\end{document}\n')



