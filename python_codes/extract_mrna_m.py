import os
import sys
import pandas as pd
import re

print("give name of gtf file for calculating gene length")
gene_len = {}
genome = open(input("enter the path of the file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'protein_coding'
    if not ((column[2] == "gene") and (re.search(pattern, str(each_line)))):
        continue
    gene = length = ''
    for elements in column[8].split(';'):
        if elements.startswith('gene_id'):
            gene = elements.replace('gene_id ', '')
            length = int(column[4])-int(column[3])
        else:
            continue
    gene_len[gene] = gene_len.get(gene, length)

print("give name of gtf file for calculating number of transcripts")
transcript = {}
genome = open(input("enter the path of the same file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'protein_coding'
    if not ((column[2] == "transcript") and (re.search(pattern, str(each_line)))):
        continue
    gene = trns = ''
    for elements in column[8].split(';'):
        if elements.startswith('gene_id'):
            gene = elements.replace('gene_id ', '')
        elif elements.startswith('transcript_id'):
            trns = elements.replace('transcript_id ', '')
        else:
            continue
    transcript[gene] = transcript.get(gene, 0) + 1

print("give name of gtf file for calculating number of exons")
exon = {}
genome = open(input("enter the path of the same file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'protein_coding'
    if not ((column[2] == "exon") and (re.search(pattern, str(each_line)))):
        continue
    gene = ex = ''
    for elements in column[8].split(';'):
        if elements.startswith('gene_id'):
            gene = elements.replace('gene_id ', '')
        elif elements.startswith('exon_id'):
            ex = elements.replace('exon_id ', '')
        else:
            continue
    exon[gene] = exon.get(gene, 0) + 1


df = pd.DataFrame([gene_len, transcript, exon])
df_t = df.T
df_t.to_csv("mrna_main_m.csv", header = False)
