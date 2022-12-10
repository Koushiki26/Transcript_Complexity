#counting number of introns for each type of splicesite

lncRNA_trns <- read.csv("lncrna_ss_m.csv", header = T)
head(lncRNA_trns)
colnames(lncRNA_trns) <- c("transcript_id", "5'SS", "3'SS", "gene_id")
length(lncRNA_trns$transcript_id)
CDS_trns <- read.table("mrna_ss_m.csv", header = T, sep = ",")
head(CDS_trns)
colnames(CDS_trns) <- c("transcript_id", "5'SS", "3'SS", "gene_id")
length(CDS_trns$transcript_id)
data_lncRNA <- read.csv("lncrna_length_m.csv", header = T)
head(data_lncRNA)
data_CDS <- read.csv("mrna_length_m.csv", header = T)
head(data_CDS)

lncRNA_trns$TPE <- ''
for(i in 1:length(lncRNA_trns$gene_id)){
  tpe <- data_lncRNA %>% filter(as.character(Gene) == as.character(lncRNA_trns$gene_id[i]))
  lncRNA_trns$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(lncRNA_trns, file = "lncrna_ss_m.csv", col.names = T, row.names = F)

CDS_trns$TPE <- ''
for(i in 1:length(CDS_trns$gene_id)){
  tpe <- data_CDS %>% filter(as.character(Gene) == as.character(CDS_trns$gene_id[i]))
  CDS_trns$TPE[i] = as.character(tpe$TPE[1])
}
write.csv(CDS_trns, file = "mrna_ss_m.csv", col.names = T, row.names = F)

lncRNA_trns <- read.csv("lncrna_ss_m.csv", header = T, sep = ",")
head(lncRNA_trns)
colnames(lncRNA_trns) <- c("transcript_id", "5'SS", "3'SS", "gene_id", "TPE")
length(lncRNA_trns$transcript_id)
lncRNA_trns <- lncRNA_trns[1:12150,]
CDS_trns <- read.table("mrna_ss_m.csv", header = T, sep = ",")
head(CDS_trns)
colnames(CDS_trns) <- c("transcript_id", "5'SS", "3'SS", "gene_id", "TPE")
length(CDS_trns$transcript_id)

x <- 0.25

a <- lncRNA_trns %>% filter(TPE <= x)
count(a)
b <- lncRNA_trns %>% filter(TPE > x)
count(b)
a1 <- CDS_trns %>% filter(TPE <= x)
count(a1)
b1 <- CDS_trns %>% filter(TPE > x)
count(b1)

#less than median 
SS <- c("aa","ac","ag","at","ca","cc","cg","ct","ga","gc","gg","gt","ta","tc","tg","tt")
y <- data.frame(SS)
y$lncRNA_count5 <- ''
y$CDS_count5 <- ''
y$lncRNA_count3 <- ''
y$CDS_count3 <- ''
for(i in 1:length(y$SS)){
  z1 <- a %>% filter(as.character(`5'SS`) == as.character(y$SS[i]))
  y$lncRNA_count5[i] = count(z1)
  z2 <- a1 %>% filter(as.character(`5'SS`) == as.character(y$SS[i]))
  y$CDS_count5[i] = count(z2)
  z3 <- a %>% filter(as.character(`3'SS`) == as.character(y$SS[i]))
  y$lncRNA_count3[i] = count(z3)
  z4 <- a1 %>% filter(as.character(`3'SS`) == as.character(y$SS[i]))
  y$CDS_count3[i] = count(z4)
}
y <- apply(y,2,as.character)
write.table(y, file = "sscount_low_TC_m.csv", col.names = T, row.names = F)

#greater than median 
SS <- c("aa","ac","ag","at","ca","cc","cg","ct","ga","gc","gg","gt","ta","tc","tg","tt")
y <- data.frame(SS)
y$lncRNA_count5 <- ''
y$CDS_count5 <- ''
y$lncRNA_count3 <- ''
y$CDS_count3 <- ''
for(i in 1:length(y$SS)){
  z1 <- b %>% filter(as.character(`5'SS`) == as.character(y$SS[i]))
  y$lncRNA_count5[i] = count(z1)
  z2 <- b1 %>% filter(as.character(`5'SS`) == as.character(y$SS[i]))
  y$CDS_count5[i] = count(z2)
  z3 <- b %>% filter(as.character(`3'SS`) == as.character(y$SS[i]))
  y$lncRNA_count3[i] = count(z3)
  z4 <- b1 %>% filter(as.character(`3'SS`) == as.character(y$SS[i]))
  y$CDS_count3[i] = count(z4)
}
y <- apply(y,2,as.character)
write.table(y, file = "sscount_high_TC_m.csv", col.names = T, row.names = F)