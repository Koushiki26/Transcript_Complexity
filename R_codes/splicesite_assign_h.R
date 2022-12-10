#assigning lncRNA and mRNA gene id for each splicesite

lncRNA_trns <- read.table("lncrna_transcript_id_h.txt", header = F, sep = "\t")
head(lncRNA_trns)
colnames(lncRNA_trns) <- c("gene_id", "transcript_id")
length(lncRNA_trns$transcript_id)
CDS_trns <- read.table("mrna_transcript_id_h.txt", header = F, sep = "\t")
head(CDS_trns)
colnames(CDS_trns) <- c("gene_id", "transcript_id")
length(CDS_trns$transcript_id)
splicesite <- read.table("splicesites_h.txt", header=F, sep="\t")
head(splicesite)
colnames(splicesite) <- c("transcript_id", "5'SS", "3'SS", "gene_id")

splicesite$gene_id <- ''
for(i in 1:length(splicesite$transcript_id)){
  ss <- lncRNA_trns %>% filter(as.character(transcript_id) == as.character(splicesite$transcript_id[i]))
  splicesite$gene_id[i] = as.character(ss$gene_id[1])
}
splicesite <- na.omit(splicesite)
write.csv(splicesite, file = "lncrna_ss_h.csv", col.names = T, row.names = F)


#run the commands again from 3-12
splicesite$gene_id <- ''
for(i in 1:length(splicesite$transcript_id)){
  ss <- CDS_trns %>% filter(as.character(transcript_id) == as.character(splicesite$transcript_id[i]))
  splicesite$gene_id[i] = as.character(ss$gene_id[1])
}
splicesite <- na.omit(splicesite)
write.csv(splicesite, file = "mrna_ss_h.csv", col.names =  T, row.names = F)
