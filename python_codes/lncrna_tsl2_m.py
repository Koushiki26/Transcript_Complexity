import os
import sys
import pandas as pd
import re

gene_len = {}
genome = open(input("Enter the path of the gene file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'lncRNA'
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

transcript = {}
genome = open(input("Enter the path of the transcript file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'lncRNA'
    if not ((column[2] == "transcript") and (re.search(pattern, str(each_line)))):
        continue
    pattern1 = 'transcript_support_level 1'
    pattern2 = 'transcript_support_level 2'
    if not ((re.search(pattern1, str(each_line))) or (re.search(pattern2, str(each_line)))):
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

exon = {}
genome = open(input("Enter the path of the exon file: "), "r")
for each_line in genome:
    if each_line.startswith('#'):
        continue
    column = each_line.strip().split('\t')
    pattern = 'lncRNA'
    if not ((column[2] == "exon") and (re.search(pattern, str(each_line)))):
        continue
    pattern1 = 'transcript_support_level 1'
    pattern2 = 'transcript_support_level 2'
    if not ((re.search(pattern1, str(each_line))) or (re.search(pattern2, str(each_line)))):
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
df_t.to_csv("lncrna_tsl2_m.csv", header = False)
