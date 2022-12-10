#installation
install.packages("devtools")
library(devtools)
install_github("monahton/GencoDymo")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
install_github("cran/seqRFLP")
remove.packages("rlang")
BiocManager::install("rlang")
devtools::install_github("r-lib/rlang", build_vignettes = TRUE)
install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.10.tar.gz", repos = NULL, type="source")
library(dplyr)
install.packages("vctrs")
install.packages("purrr")

#maxentscan tool analysis
library(GencoDymo)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome)
df <- load_gtf("gencode.v38.basic.annotation.gtf")
head(df)
intron_df <- extract_introns(df)
head(intron_df)
library(dplyr)
data_lncRNA <- read.csv("lncRNA_length.csv", header = T)
head(data_lncRNA)
data_CDS <- read.csv("CDS_length.csv", header = T)
head(data_CDS)

#5ss
#for lncRNA
introns_df <- intron_df %>% filter(transcript_type=="lncRNA") %>% distinct(intron_start, intron_end, .keep_all = T)
head(introns_df)
introns_df$TPE <- ''
for(i in 1:length(introns_df$gene_id)){
  tpe <- data_lncRNA %>% filter(as.character(Gene) == as.character(introns_df$gene_id[i]))
  introns_df$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(introns_df, file = "lncRNA_max_h.csv", col.names = T, row.names = F)
x <- 0.3333333
a <- introns_df %>% filter(TPE <= x)
count(a)
extract_5ss_motif(a, genome = BSgenome.Hsapiens.UCSC.hg38)
b <- introns_df %>% filter(TPE > x)
count(b)
extract_5ss_motif(b, genome = BSgenome.Hsapiens.UCSC.hg38)

#for mRNA
introns_df <- intron_df %>% filter(transcript_type=="protein_coding") %>% distinct(intron_start, intron_end, .keep_all = T)
head(introns_df)
introns_df$TPE <- ''
for(i in 1:length(introns_df$gene_id)){
  tpe <- data_CDS %>% filter(as.character(Gene) == as.character(introns_df$gene_id[i]))
  introns_df$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(introns_df, file = "CDS_max_h.csv", col.names = T, row.names = F)
x <- 0.3333333
a <- introns_df %>% filter(TPE <= x)
count(a)
extract_5ss_motif(a, genome = BSgenome.Hsapiens.UCSC.hg38)
b <- introns_df %>% filter(TPE > x)
count(b)
extract_5ss_motif(b, genome = BSgenome.Hsapiens.UCSC.hg38)


#3ss
#for lncRNA
introns_df <- intron_df %>% filter(transcript_type=="lncRNA") %>% distinct(intron_start, intron_end, .keep_all = T)
head(introns_df)
introns_df$TPE <- ''
for(i in 1:length(introns_df$gene_id)){
  tpe <- data_lncRNA %>% filter(as.character(Gene) == as.character(introns_df$gene_id[i]))
  introns_df$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(introns_df, file = "lncRNA_max_h.csv", col.names = T, row.names = F)
x <- 0.3333333
a <- introns_df %>% filter(TPE <= x)
count(a)
extract_3ss_motif(a, genome = BSgenome.Hsapiens.UCSC.hg38)
b <- introns_df %>% filter(TPE > x)
count(b)
extract_3ss_motif(b, genome = BSgenome.Hsapiens.UCSC.hg38)

#for mRNA
introns_df <- intron_df %>% filter(transcript_type=="protein_coding") %>% distinct(intron_start, intron_end, .keep_all = T)
head(introns_df)
introns_df$TPE <- ''
for(i in 1:length(introns_df$gene_id)){
  tpe <- data_CDS %>% filter(as.character(Gene) == as.character(introns_df$gene_id[i]))
  introns_df$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(introns_df, file = "CDS_max_h.csv", col.names = T, row.names = F)
x <- 0.3333333
a <- introns_df %>% filter(TPE <= x)
count(a)
extract_3ss_motif(a, genome = BSgenome.Hsapiens.UCSC.hg38)
b <- introns_df %>% filter(TPE > x)
count(b)
extract_3ss_motif(b, genome = BSgenome.Hsapiens.UCSC.hg38)

