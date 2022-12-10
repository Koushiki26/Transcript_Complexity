#assigning the intron and and exon length against lncRNA and mRNA gene id

library(dplyr)
data_intron <- read.table("intron_length_h.txt", header = F, sep = "\t")
head(data_intron)
colnames(data_intron) <- c("Transcript_id", "intron_length")
data_exon <- read.table("exon_length_h.txt", header = F, sep = "\t")
head(data_exon)
colnames(data_exon) <- c("Transcript_id", "exon_length")
lncRNA_trns <- read.table("lncrna_transcript_id_h.txt", header = F, sep = "\t")
head(lncRNA_trns)
colnames(lncRNA_trns) <- c("gene_id", "transcript_id")
length(lncRNA_trns$transcript_id)
CDS_trns <- read.table("mrna_transcript_id_h.txt", header = F, sep = "\t")
head(CDS_trns)
colnames(CDS_trns) <- c("gene_id", "transcript_id")
length(CDS_trns$transcript_id)

lncRNA_trns$intron_length <- ''
lncRNA_trns$exon_length <- ''
for(i in 1:length(lncRNA_trns$transcript_id))
{
  intron <- data_intron %>% filter(as.character(Transcript_id) == as.character(lncRNA_trns$transcript_id[i]))
  lncRNA_trns$intron_length[i] = intron$intron_length[1]
  exon <- data_exon %>% filter(as.character(Transcript_id) == as.character(lncRNA_trns$transcript_id[i]))
  lncRNA_trns$exon_length[i] = exon$exon_length[1]
}
summary(lncRNA_trns)

CDS_trns$intron_length <- ''
CDS_trns$exon_length <- ''
for(i in 1:length(CDS_trns$transcript_id))
{
  intron <- data_intron %>% filter(as.character(Transcript_id) == as.character(CDS_trns$transcript_id[i]))
  CDS_trns$intron_length[i] = intron$intron_length[1]
  exon <- data_exon %>% filter(as.character(Transcript_id) == as.character(CDS_trns$transcript_id[i]))
  CDS_trns$exon_length[i] = exon$exon_length[1]
}
summary(CDS_trns)

write.csv(lncRNA_trns, file = "lncrna_trns_h.csv", col.names = T, row.names = F, na = "0")
write.csv(CDS_trns, file = "mrna_trns_h.csv", col.names = T, row.names = F, na = "0")


data_lncRNA <- read.csv("lncrna_main_h.csv", header = F)
colnames(data_lncRNA) <- c('Gene','Gene_length','Transcript','Exon')
head(data_lncRNA)
data_CDS <- read.csv("mrna_main_h.csv", header = F)
colnames(data_CDS) <- c('Gene','Gene_length','Transcript','Exon')
head(data_CDS)
lncRNA_trns <- read.csv("lncrna_trns_h.csv", header = T)
head(lncRNA_trns)
length(lncRNA_trns$transcript_id)
CDS_trns <- read.csv("mrna_trns_h.csv", header = T)
head(CDS_trns)
length(CDS_trns$transcript_id)

summary(lncRNA_trns)
lncRNA_trns$intron_length <- as.numeric(lncRNA_trns$intron_length)
lncRNA_trns$exon_length <- as.numeric(lncRNA_trns$exon_length)
summary(CDS_trns)
CDS_trns$intron_length <- as.numeric(CDS_trns$intron_length)
CDS_trns$exon_length <- as.numeric(CDS_trns$exon_length)

data_lncRNA$intron_length <- ''
data_lncRNA$exon_length <- ''
for(i in 1:length(data_lncRNA$Gene))
{
  intron <- lncRNA_trns %>% filter(as.character(gene_id) == as.character(data_lncRNA$Gene[i]))
  data_lncRNA$intron_length[i] = sum(intron$intron_length)
  exon <- lncRNA_trns %>% filter(as.character(gene_id) == as.character(data_lncRNA$Gene[i]))
  data_lncRNA$exon_length[i] = sum(exon$exon_length)
}

data_CDS$intron_length <- ''
data_CDS$exon_length <- ''
for(i in 1:length(data_CDS$Gene))
{
  intron <- CDS_trns %>% filter(as.character(gene_id) == as.character(data_CDS$Gene[i]))
  data_CDS$intron_length[i] = sum(intron$intron_length)
  exon <- CDS_trns %>% filter(as.character(gene_id) == as.character(data_CDS$Gene[i]))
  data_CDS$exon_length[i] = sum(exon$exon_length)
}

data_lncRNA$intron_length <- as.numeric(data_lncRNA$intron_length)
data_lncRNA$exon_length <- as.numeric(data_lncRNA$exon_length)
summary(data_lncRNA)
data_CDS$intron_length <- as.numeric(data_CDS$intron_length)
data_CDS$exon_length <- as.numeric(data_CDS$exon_length)
summary(data_CDS)
data_lncRNA <- mutate(data_lncRNA, TPE = Transcript/Exon)
data_CDS <- mutate(data_CDS, TPE = Transcript/Exon)

write.csv(data_lncRNA, file = "lncrna_length_h.csv", col.names = T, row.names = F, na = "0")
write.csv(data_CDS, file = "mrna_length_h.csv", col.names = T, row.names = F, na = "0")