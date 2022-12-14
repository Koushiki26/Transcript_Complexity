# Transcript_Complexity
Genome-wide analysis of lncRNAs and mRNAs Transcript Complexity in human and mouse

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
1. as_h_A3_strict.ioe :
2. as_h_A5_strict.ioe :
3. as_h_AF_strict.ioe :
4. as_h_AL_strict.ioe :
5. as_h_MX_strict.ioe :
6. as_h_RI_strict.ioe :
7. as_h_SE_strict.ioe :
8. as_m_A3_strict.ioe :
9. as_m_A5_strict.ioe :
10. as_m_AF_strict.ioe :
11. as_m_AL_strict.ioe :
12. as_m_MX_strict.ioe :
13. as_m_RI_strict.ioe :
14. as_m_SE_strict.ioe :
15. eclip_annotation_filtered_h.txt :
16. exon_length_h.txt :
17. exon_length_m.txt :
18. final_lncrna_high_3sss_h.txt :
19. final_lncrna_high_5sss_h.txt :
20. final_lncrna_low_3sss_h.txt :
21. final_lncrna_low_5sss_h.txt :
22. final_mrna_high_3sss_h.txt :
23. final_mrna_high_5sss_h.txt :
24. final_mrna_low_3sss_h.txt :
25. final_mrna_low_5sss_h.txt :
26. intron_length_h.txt :
27. intron_length_m.txt :
28. lncrna_ASE_h.csv :
29. lncrna_ASE_m.csv :
30. lncrna_length_h.csv :
31. lncrna_length_m.csv :
32. lncrna_main_h.csv :
33. lncrna_main_m.csv :
34. lncrna_ss_h.csv :
35. lncrna_ss_m.csv :
36. lncrna_transcript_id_h.txt :
37. lncrna_transcript_id_m.txt :
38. lncrna_trns_h.csv :
39. lncrna_trns_m.csv :
40. lncrna_tsl1_h.csv :
41. lncrna_tsl1_m.csv :
42. lncrna_tsl2_h.csv :
43. lncrna_tsl2_m.csv :
44. lncrna_tsl3_h.csv :
45. lncrna_tsl3_m.csv :
46. lncrna_tsl4_h.csv :
47. lncrna_tsl4_m.csv :
48. lncrna_tsl5_h.csv :
49. lncrna_tsl5_m.csv :
50. lncrna_tsl6_h.csv :
51. lncrna_tsl6_m.csv :
52. mrna_ASE_h.csv :
53. mrna_ASE_m.csv :
54. mrna_eclip_h.csv :
55. mrna_length_h.csv :
56. mrna_length_m.csv :
57. mrna_main_h.csv :
58. mrna_main_m.csv :
59. mrna_ss_h.csv :
60. mrna_ss_m.csv :
61. mrna_transcript_id_h.txt :
62. mrna_transcript_id_m.txt :
63. mrna_trns_h.csv :
64. mrna_trns_m.csv :
65. mrna_tsl1_h.csv :
66. mrna_tsl1_m.csv :
67. mrna_tsl2_h.csv :
68. mrna_tsl2_m.csv :
69. mrna_tsl3_h.csv :
70. mrna_tsl3_m.csv :
71. mrna_tsl4_h.csv :
72. mrna_tsl4_m.csv :
73. mrna_tsl5_h.csv :
74. mrna_tsl5_m.csv :
75. mrna_tsl6_h.csv :
76. mrna_tsl6_m.csv :
77. splicesites_h.txt :
78. splicesites_m.txt :
79. sscount_high_TC_h.csv :
80. sscount_high_TC_m.csv :
81. sscount_low_TC_h.csv :
82. sscount_low_TC_m.csv :

--------------------------------
# data_human_conservation
--------------------------------

--------------------------------
# data_mouse_conservation
--------------------------------

--------------------------------
# supplementary_files
--------------------------------
1. SUPPLEMENTARY MATERIAL.pdf :
2. supplemental_file1.xlsx :
3. supplemental_file2.xlsx :
4. supplemental_file3.xlsx :
5. supplemental_file4.xlsx :
