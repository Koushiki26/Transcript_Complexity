# Transcript_Complexity
Genome-wide analysis of lncRNAs and mRNAs Transcript Complexity in human and mouse

----------------------------
# Table of Contents
----------------------------

   * [Python codes](#python_codes)
   * [R codes](#R_codes)
   * [Data](#data)
   * [Human Conservation Data](#data_human_conservation)
   * [Mouse Conservation Data](#data_mouse_conservation)
   * [Supplementary files](#supplementary_files)

--------------------------------
# python_codes
--------------------------------
1. coordinate.py : Take FASTA file of intron sequences as input and extract the IDs (containing coordinates) of each sequence.
2. extract_lncrna_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene.
3. extract_lncrna_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene.
4. extract_mrna_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene.
5. extract_mrna_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene.
6. length.py : Take fasta file of exon and intron sequences as input and calculate the length for each sequence.
7. lncrna_tsl1_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1.
8. lncrna_tsl1_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1.
9. lncrna_tsl2_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1 or 2.
10. lncrna_tsl2_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1 or 2.
11. lncrna_tsl3_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2 or 3.
12. lncrna_tsl3_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2 or 3.
13. lncrna_tsl4_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3 or 4.
14. lncrna_tsl4_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3 or 4.
15. lncrna_tsl5_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3, 4 or 5.
16. lncrna_tsl5_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3, 4 or 5.
17. lncrna_tsl6_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3, 4, 5 or NA.
18. lncrna_tsl6_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each lncRNA gene with Transcript support level 1, 2, 3, 4, 5 or NA.
19. mrna_tsl1_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1.
20. mrna_tsl1_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1.
21. mrna_tsl2_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1 or 2.
22. mrna_tsl2_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1 or 2.
23. mrna_tsl3_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2 or 3.
24. mrna_tsl3_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2 or 3.
25. mrna_tsl4_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3 or 4.
26. mrna_tsl4_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3 or 4.
27. mrna_tsl5_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3, 4 or 5.
28. mrna_tsl5_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3, 4 or 5.
29. mrna_tsl6_h.py : Take gtf file (human) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3, 4, 5 or NA.
30. mrna_tsl6_m.py : Take gtf file (mouse) as input and calculate the gene length, number of transcript and number of exon for each mRNA gene with Transcript support level 1, 2, 3, 4, 5 or NA.
31. splicesite.py : Takes FASTA file of intron sequences as input and extract the 5' and 3' splice sites of each sequence.

--------------------------------
# R_codes
--------------------------------
1. ASE_h.R : Take ioe files (human) as input and calculate the number of transcripts involved in different alternative splicing event for each lncRNA and mRNA gene.
2. ASE_m.R : Take ioe files (mouse) as input and calculate the number of transcripts involved in different alternative splicing event for each lncRNA and mRNA gene.
3. ASE_plots.R : Take the alternative splicing event data for lncRNA and mRNA genes as input and plot the fraction of transcripts involved in different alternative splicing event (human and mouse).
4. conservation_score_plot_h.R : Take lncRNA (low and high TC) and mRNA (low and high TC) conservation score data as input and plot the graph for 5' and 3' splice site conservation (human).
5. conservation_score_plot_m.R : Take lncRNA (low and high TC) and mRNA (low and high TC) conservation score data as input and plot the graph for 5' and 3' splice site conservation (mouse).
6. correlation_analysis_h.R : Take the number of transcript and exon data of lncRNA and mRNA genes as the input and calculate and plot correlation between them (human).
7. correlation_analysis_m.R : Take the number of transcript and exon data of lncRNA and mRNA genes as the input and calculate and plot correlation between them (mouse).
8. eclip.R : Take the number of transcript and exon data of lncRNA and mRNA genes (NMD analysis) as the input and calculate and plot correlation between them (human).
9. length_h.R : Assign the intron and exon length against lncRNA and mRNA gene ID (human).
10. length_m.R : Assign the intron and exon length against lncRNA and mRNA gene ID (mouse).
11. length_plot_h.R : Take exon and inron length data of lncRNA and mRNA genes as input and plot length data with respect to transcript complexity (human).
12. length_plot_m.R : Take exon and inron length data of lncRNA and mRNA genes as input and plot length data with respect to transcript complexity (mouse).
13. maxentscan_h.R : Create FASTA files of lncRNA and mRNA genes with high and low Transcript Complexity for input in MaxEntScan tool (human).
14. mean_TC.R : Take the number of transcript and exon data of lncRNA and mRNA genes as the input and calculate and plot Transcript Complexity (human and mouse).
15. mean_TC_tsl_h.R : Take the number of transcript and exon data of lncRNA and mRNA genes for different Transcript Support Level as the input and calculate and plot Transcript Complexity (human).
16. mean_TC_tsl_m.R : Take the number of transcript and exon data of lncRNA and mRNA genes for different Transcript Support Level as the input and calculate and plot Transcript Complexity (mouse).
17. splicesite_assign_h.R : Take the 5' and 3' splicesite dinucleotide data for each intron as input and assign the lncRNA and mRNA gene ID to each splicesite (human).
18. splicesite_assign_m.R : Take the 5' and 3' splicesite dinucleotide data for each intron as input and assign the lncRNA and mRNA gene ID to each splicesite (mouse).
19. splicesite_count_h.R : Count number of introns for each type of splicesite dinucleotide and then sort them based on high or low Transcript Complexity (human).
20. splicesite_count_m.R : Count number of introns for each type of splicesite dinucleotide and then sort them based on high or low Transcript Complexity (mouse).
21. splicesite_plot_h.R : Take the splicesite count data as input and plot the fraction of intron for each 5' and 3' splicesite dinucleotide (human).
22. splicesite_plot_m.R : Take the splicesite count data as input and plot the fraction of intron for each 5' and 3' splicesite dinucleotide (mouse).
23. splicesite_strength_plot_h.R : Take the splicesite strength data as input and plot it (human).

--------------------------------
# data
--------------------------------
1. as_h_A3_strict.ioe : Alternative 3' splice site data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
2. as_h_A5_strict.ioe : Alternative 5' splice site data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
3. as_h_AF_strict.ioe : Alternative First Exons data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
4. as_h_AL_strict.ioe : Alternative Last Exons data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
5. as_h_MX_strict.ioe : Mutually Exclusive Exons data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
6. as_h_RI_strict.ioe : Retained Intron data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
7. as_h_SE_strict.ioe : Skipping Exons data of human (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
8. as_m_A3_strict.ioe : Alternative 3' splice site data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
9. as_m_A5_strict.ioe : Alternative 5' splice site data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
10. as_m_AF_strict.ioe : Alternative First Exons data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
11. as_m_AL_strict.ioe : Alternative Last Exons data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
12. as_m_MX_strict.ioe : Mutually Exclusive Exons data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
13. as_m_RI_strict.ioe : Retained Intron data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
14. as_m_SE_strict.ioe : Skipping Exons data of mouse (column names-seqname, gene id, event id, alternative transcripts, total transcripts)
15. eclip_annotation_filtered_h.txt : Annotated eCLIP data of human obtained by using annotatePeaks.pl command of HOMER (column names-peak id, chr, start position, stop position, strand, gene id, type of gene)
16. exon_length_h.txt : total exon length data for each transcript of human (column names-transcript id, exon length)
17. exon_length_m.txt : total exon length data for each transcript of mouse (column names-transcript id, exon length)
18. final_lncrna_high_3sss_h.txt : MaxEntScan::score3ss or 3' splice site strength data for lncRNA with high TC of human (column names-sequence, MAXENT, MM, WMM)
19. final_lncrna_high_5sss_h.txt : MaxEntScan::score5ss or 5' splice site strength data for lncRNA with high TC of human (column names-sequence, MAXENT, MDD, MM, WMM)
20. final_lncrna_low_3sss_h.txt : MaxEntScan::score3ss or 3' splice site strength data for lncRNA with low TC of human (column names-sequence, MAXENT, MM, WMM)
21. final_lncrna_low_5sss_h.txt : MaxEntScan::score5ss or 5' splice site strength data for lncRNA with low TC of human (column names-sequence, MAXENT, MDD, MM, WMM)
22. final_mrna_high_3sss_h.txt : MaxEntScan::score3ss or 3' splice site strength data for mRNA with high TC of human (column names-sequence, MAXENT, MM, WMM)
23. final_mrna_high_5sss_h.txt : MaxEntScan::score5ss or 5' splice site strength data for mRNA with high TC of human (column names-sequence, MAXENT, MDD, MM, WMM)
24. final_mrna_low_3sss_h.txt : MaxEntScan::score3ss or 3' splice site strength data for mRNA with low TC of human (column names-sequence, MAXENT, MM, WMM)
25. final_mrna_low_5sss_h.txt : MaxEntScan::score5ss or 5' splice site strength data for mRNA with low TC of human (column names-sequence, MAXENT, MDD, MM, WMM)
26. intron_length_h.txt : total intron length data for each transcript of human (column names-transcript id, intron length)
27. intron_length_m.txt : total intron length data for each transcript of mouse (column names-transcript id, intron length)
28. lncrna_ASE_h.csv : count of transcripts data in different Alternative Splicing Events for each lncRNA genes of human (column names-gene id, gene length, number of transcripts, number of exons, A3, A5, AF, AL, MX, RI, SE)
29. lncrna_ASE_m.csv : count of transcripts data in different Alternative Splicing Events for each lncRNA genes of mouse (column names-gene id, gene length, number of transcripts, number of exons, A3, A5, AF, AL, MX, RI, SE)
30. lncrna_length_h.csv : exon length and intron length data for each lncRNA gene of human (column names-gene id, gene length, number of transcripts, number of exons, intron length, exon length, TC)
31. lncrna_length_m.csv : exon length and intron length data for each lncRNA gene of mouse (column names-gene id, gene length, number of transcripts, number of exons, intron length, exon length, TC)
32. lncrna_main_h.csv : number of transcript and number of exon data for each lncRNA gene of human (column names-gene id, gene length, number of transcripts, number of exons)
33. lncrna_main_m.csv : number of transcript and number of exon data for each lncRNA gene of mouse (column names-gene id, gene length, number of transcripts, number of exons)
34. lncrna_ss_h.csv : 5' and 3' splicesite data for each lncRNA transcript of human (column names-transcript id, 5' splicesite, 3' splicesite, gene id, TC)
35. lncrna_ss_m.csv : 5' and 3' splicesite data for each lncRNA transcript of mouse (column names-transcript id, 5' splicesite, 3' splicesite, gene id, TC)
36. lncrna_transcript_id_h.txt : transcript id data for lncRNA genes of human (column names-gene id, transcript id)
37. lncrna_transcript_id_m.txt : transcript id data for lncRNA genes of mouse (column names-gene id, transcript id)
38. lncrna_trns_h.csv : exon length and intron length data for each lncRNA transcript of human (column names-gene id, transcript id, intron length, exon length)
39. lncrna_trns_m.csv : exon length and intron length data for each lncRNA transcript of mouse (column names-gene id, transcript id, intron length, exon length)
40. lncrna_tsl1_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1 of human (column names-gene id, gene length, number of transcripts, number of exons)
41. lncrna_tsl1_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
42. lncrna_tsl2_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1 or 2 of human (column names-gene id, gene length, number of transcripts, number of exons)
43. lncrna_tsl2_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1 or 2 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
44. lncrna_tsl3_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2 or 3 of human (column names-gene id, gene length, number of transcripts, number of exons)
45. lncrna_tsl3_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2 or 3 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
46. lncrna_tsl4_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3 or 4 of human (column names-gene id, gene length, number of transcripts, number of exons)
47. lncrna_tsl4_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3 or 4 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
48. lncrna_tsl5_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3, 4 or 5 of human (column names-gene id, gene length, number of transcripts, number of exons)
49. lncrna_tsl5_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3, 4 or 5 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
50. lncrna_tsl6_h.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3, 4, 5 or NA of human (column names-gene id, gene length, number of transcripts, number of exons)
51. lncrna_tsl6_m.csv : number of transcript and number of exon data for each lncRNA gene with Transcript Support Level 1, 2, 3, 4, 5 or NA of mouse (column names-gene id, gene length, number of transcripts, number of exons)
52. mrna_ASE_h.csv : count of transcripts data in different Alternative Splicing Events for each mRNA genes of human (column names-gene id, gene length, number of transcripts, number of exons, A3, A5, AF, AL, MX, RI, SE)
53. mrna_ASE_m.csv : count of transcripts data in different Alternative Splicing Events for each mRNA genes of mouse (column names-gene id, gene length, number of transcripts, number of exons, A3, A5, AF, AL, MX, RI, SE)
54. mrna_eclip_h.csv : number of transcript and number of exon data for each mRNA gene from eCLIP data of human (column names-gene id, gene length, number of transcripts, number of exons)
55. mrna_length_h.csv : exon length and intron length data for each mRNA gene of human (column names-gene id, gene length, number of transcripts, number of exons, intron length, exon length, TC)
56. mrna_length_m.csv : exon length and intron length data for each mRNA gene of mouse (column names-gene id, gene length, number of transcripts, number of exons, intron length, exon length, TC)
57. mrna_main_h.csv : number of transcript and number of exon data for each mRNA gene of human (column names-gene id, gene length, number of transcripts, number of exons)
58. mrna_main_m.csv : number of transcript and number of exon data for each mRNA gene of mouse (column names-gene id, gene length, number of transcripts, number of exons)
59. mrna_ss_h.csv : 5' and 3' splicesite data for each mRNA transcript of human (column names-transcript id, 5' splicesite, 3' splicesite, gene id, TC)
60. mrna_ss_m.csv : 5' and 3' splicesite data for each mRNA transcript of mouse (column names-transcript id, 5' splicesite, 3' splicesite, gene id, TC)
61. mrna_transcript_id_h.txt : transcript id data for mRNA genes of human (column names-gene id, transcript id)
62. mrna_transcript_id_m.txt : transcript id data for mRNA genes of mouse (column names-gene id, transcript id)
63. mrna_trns_h.csv : exon length and intron length data for each mRNA transcript of human (column names-gene id, transcript id, intron length, exon length)
64. mrna_trns_m.csv : exon length and intron length data for each mRNA transcript of mouse (column names-gene id, transcript id, intron length, exon length)
65. mrna_tsl1_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1 of human (column names-gene id, gene length, number of transcripts, number of exons)
66. mrna_tsl1_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
67. mrna_tsl2_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1 or 2 of human (column names-gene id, gene length, number of transcripts, number of exons)
68. mrna_tsl2_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1 or 2 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
69. mrna_tsl3_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2 or 3 of human (column names-gene id, gene length, number of transcripts, number of exons)
70. mrna_tsl3_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2 or 3 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
71. mrna_tsl4_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3 or 4 of human (column names-gene id, gene length, number of transcripts, number of exons)
72. mrna_tsl4_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3 or 4 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
73. mrna_tsl5_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3, 4 or 5 of human (column names-gene id, gene length, number of transcripts, number of exons)
74. mrna_tsl5_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3, 4 or 5 of mouse (column names-gene id, gene length, number of transcripts, number of exons)
75. mrna_tsl6_h.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3, 4, 5 or NA of human (column names-gene id, gene length, number of transcripts, number of exons)
76. mrna_tsl6_m.csv : number of transcript and number of exon data for each mRNA gene with Transcript Support Level 1, 2, 3, 4, 5 or NA of mouse (column names-gene id, gene length, number of transcripts, number of exons)
77. splicesites_h.txt : 5' and 3' splicesite for each intron of human (column names-transcript id, 5' splicesite, 3' splicesite)
78. splicesites_m.txt : 5' and 3' splicesite for each intron of mouse (column names-transcript id, 5' splicesite, 3' splicesite)
79. sscount_high_TC_h.csv : count of each 5' and 3' splicesite dinucleotide for lncRNA and mRNA genes with high Transcript complexity of human (column names-splicesite, 5' lncRNA count, 5' mRNA count,  3' lncRNA count, 3' mRNA count)
80. sscount_high_TC_m.csv : count of each 5' and 3' splicesite dinucleotide for lncRNA and mRNA genes with high Transcript complexity of mouse (column names-splicesite, 5' lncRNA count, 5' mRNA count,  3' lncRNA count, 3' mRNA count)
81. sscount_low_TC_h.csv : count of each 5' and 3' splicesite dinucleotide for lncRNA and mRNA genes with low Transcript complexity of human (column names-splicesite, 5' lncRNA count, 5' mRNA count,  3' lncRNA count, 3' mRNA count)
82. sscount_low_TC_m.csv : count of each 5' and 3' splicesite dinucleotide for lncRNA and mRNA genes with low Transcript complexity of mouse (column names-splicesite, 5' lncRNA count, 5' mRNA count,  3' lncRNA count, 3' mRNA count)

--------------------------------
# data_human_conservation
--------------------------------
This folder contains the conservation score of human 5' and 3' splicesite for lncRNA and mRNA genes with high and low Transcript complexity for every chromosome. column names - chromosome number, position, phastcon score

--------------------------------
# data_mouse_conservation
--------------------------------
This folder contains the conservation score of mouse 5' and 3' splicesite for lncRNA and mRNA genes with high and low Transcript complexity for every chromosome. column names - chromosome number, position, phastcon score

--------------------------------
# supplementary_files
--------------------------------
1. SUPPLEMENTARY MATERIAL.pdf : contains all the supplementary tables and supplementary figures
2. supplemental_file1.xlsx : contain data of human lncRNA and mRNA genes
a) lncRNA gene id, gene length, number of transcript, number of exon, transcript complexity
b) lncRNA gene id, A3, A5, AF, AL, MX, RI, SE
c) lncRNA gene id, intron length, exon length
d) lncRNA transcript id, 5' splicesite, 3' splicesite, lncRNA gene id
e) mRNA gene id, gene length, number of transcript, number of exon, transcript complexity
f) mRNA gene id, A3, A5, AF, AL, MX, RI, SE
g) mRNA gene id, intron length, exon length
h) mRNA transcript id, 5' splicesite, 3' splicesite, mRNA gene id
3. supplemental_file2.xlsx : contain data of mouse lncRNA and mRNA genes
a) lncRNA gene id, gene length, number of transcript, number of exon, transcript complexity
b) lncRNA gene id, A3, A5, AF, AL, MX, RI, SE
c) lncRNA gene id, intron length, exon length
d) lncRNA transcript id, 5' splicesite, 3' splicesite, lncRNA gene id
e) mRNA gene id, gene length, number of transcript, number of exon, transcript complexity
f) mRNA gene id, A3, A5, AF, AL, MX, RI, SE
g) mRNA gene id, intron length, exon length
h) mRNA transcript id, 5' splicesite, 3' splicesite, mRNA gene id
4. supplemental_file3.xlsx : number of transcript and number of exon data of human with different Transcript Support Level
a) lncRNA (TSL1)
b) lncRNA (TSL1+TSL2)
c) lncRNA (TSL1+TSL2+TSL3)
d) lncRNA (TSL1+TSL2+TSL3+TSL4)
e) lncRNA (TSL1+TSL2+TSL3+TSL4+TSL5)
f) lncRNA (TSL1+TSL2+TSL3+TSL4+TSL5+TSLNA)
g) mRNA (TSL1)
h) mRNA (TSL1+TSL2)
i) mRNA (TSL1+TSL2+TSL3)
j) mRNA (TSL1+TSL2+TSL3+TSL4)
k) mRNA (TSL1+TSL2+TSL3+TSL4+TSL5)
l) mRNA (TSL1+TSL2+TSL3+TSL4+TSL5+TSLNA)
6. supplemental_file4.xlsx : number of transcript and number of exon data of mouse with different Transcript Support Level
a) lncRNA (TSL1)
b) lncRNA (TSL1+TSL2)
c) lncRNA (TSL1+TSL2+TSL3)
d) lncRNA (TSL1+TSL2+TSL3+TSL4)
e) lncRNA (TSL1+TSL2+TSL3+TSL4+TSL5)
f) lncRNA (TSL1+TSL2+TSL3+TSL4+TSL5+TSLNA)
g) mRNA (TSL1)
h) mRNA (TSL1+TSL2)
i) mRNA (TSL1+TSL2+TSL3)
j) mRNA (TSL1+TSL2+TSL3+TSL4)
k) mRNA (TSL1+TSL2+TSL3+TSL4+TSL5)
l) mRNA (TSL1+TSL2+TSL3+TSL4+TSL5+TSLNA)
