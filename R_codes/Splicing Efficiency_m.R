#splicing efficiency
data_lncRNA <- read.csv("lncrna_main_m.csv", header = F)
head(data_lncRNA)
colnames(data_lncRNA) <- c("gene_ID", "gene_length", "transcript","exon")
data_lncRNA$TC <- data_lncRNA$transcript/data_lncRNA$exon
data_CDS <- read.csv("mrna_main_m.csv", header = F)
head(data_CDS)
colnames(data_CDS) <- c("gene_ID", "gene_length", "transcript","exon")
data_CDS$TC <- data_CDS$transcript/data_CDS$exon
G1E_data <- read.table("G1E_output1.tsv", header = F, sep="\t")
colnames(G1E_data) <- c("chr", "strand", "gene_ID", "transcript_ID","intron_ID", 
                        "sj5start","sj5end","sj5_cov_nonsplit","sj5_cov_split",
                        "sj3start", "sj3end", "sj3_cov_nonsplit","sj3_cov_split","score")
head(G1E_data)
summary(G1E_data)
MEL_data <- read.table("MEL_output1.tsv", header = F, sep = "\t")
colnames(MEL_data) <- c("chr", "strand", "gene_ID", "transcript_ID","intron_ID", 
                        "sj5start","sj5end","sj5_cov_nonsplit","sj5_cov_split",
                        "sj3start", "sj3end", "sj3_cov_nonsplit","sj3_cov_split","score")
head(MEL_data)
summary(MEL_data)

library(dplyr)
G1E_data$TC <- ''
G1E_data$exon <- ''
for(i in 1:length(G1E_data$score))
{
  TC <- data_lncRNA %>% filter(as.character(gene_ID) == as.character(G1E_data$gene_ID[i]))
  G1E_data$TC[i] = TC$TC[1]
  G1E_data$exon[i] = TC$exon[1]
}
G1E_data <- na.omit(G1E_data)
write.csv(G1E_data, file = "lncRNA_G1E.csv", col.names = T, row.names = F)

G1E_data$TC <- ''
G1E_data$exon <- ''
for(i in 1:length(G1E_data$score))
{
  TC <- data_CDS %>% filter(as.character(gene_ID) == as.character(G1E_data$gene_ID[i]))
  G1E_data$TC[i] = TC$TC[1]
  G1E_data$exon[i] = TC$exon[1]
}
G1E_data <- na.omit(G1E_data)
write.csv(G1E_data, file = "mRNA_G1E.csv", col.names = T, row.names = F)

data_lncRNA <- read.csv("lncRNA_G1E.csv", header = T)
head(data_lncRNA)
data_pcg <- read.csv("mRNA_G1E.csv", header = T)
head(data_pcg)

library(dplyr)
a <- data_lncRNA %>% filter(exon<11)
count(a)
b <- data_lncRNA %>% filter(exon>10&exon<21)
count(b)
c <- data_lncRNA %>% filter(exon>20&exon<51)
count(c)
d <- data_lncRNA %>% filter(exon>50&exon<101)
count(d)
e <- data_lncRNA %>% filter(exon>100)
count(e)
a1 <- data_pcg %>% filter(exon<11)
count(a1)
b1 <- data_pcg %>% filter(exon>10&exon<21)
count(b1)
c1 <- data_pcg %>% filter(exon>20&exon<51)
count(c1)
d1 <- data_pcg %>% filter(exon>50&exon<101)
count(d1)
e1 <- data_pcg %>% filter(exon>100)
count(e1)
tiff("/home/dell/Documents/koushiki/final_fig/SE_bin1_G1E.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$score, b$score, c$score, d$score, e$score,
        col = "#4D648D", xlab="Number of Exon", ylab="Splicing Efficiency", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), main="")
boxplot(a1$score, b1$score, c1$score, d1$score, e1$score,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()

MEL_data$TC <- ''
MEL_data$exon <- ''
for(i in 1:length(MEL_data$score))
{
  TC <- data_lncRNA %>% filter(as.character(gene_ID) == as.character(MEL_data$gene_ID[i]))
  MEL_data$TC[i] = TC$TC[1]
  MEL_data$exon[i] = TC$exon[1]
}
MEL_data <- na.omit(MEL_data)
write.csv(MEL_data, file = "lncRNA_MEL.csv", col.names = T, row.names = F)

MEL_data$TC <- ''
MEL_data$exon <- ''
for(i in 1:length(MEL_data$score))
{
  TC <- data_CDS %>% filter(as.character(gene_ID) == as.character(MEL_data$gene_ID[i]))
  MEL_data$TC[i] = TC$TC[1]
  MEL_data$exon[i] = TC$exon[1]
}
MEL_data <- na.omit(MEL_data)
write.csv(MEL_data, file = "mRNA_MEL.csv", col.names = T, row.names = F)

data_lncRNA <- read.csv("lncRNA_MEL.csv", header = T)
head(data_lncRNA)
data_pcg <- read.csv("mRNA_MEL.csv", header = T)
head(data_pcg)

a <- data_lncRNA %>% filter(exon<11)
count(a)
b <- data_lncRNA %>% filter(exon>10&exon<21)
count(b)
c <- data_lncRNA %>% filter(exon>20&exon<51)
count(c)
d <- data_lncRNA %>% filter(exon>50&exon<101)
count(d)
e <- data_lncRNA %>% filter(exon>100)
count(e)
a1 <- data_pcg %>% filter(exon<11)
count(a1)
b1 <- data_pcg %>% filter(exon>10&exon<21)
count(b1)
c1 <- data_pcg %>% filter(exon>20&exon<51)
count(c1)
d1 <- data_pcg %>% filter(exon>50&exon<101)
count(d1)
e1 <- data_pcg %>% filter(exon>100)
count(e1)
tiff("/home/dell/Documents/koushiki/final_fig/SE_bin1_MEL.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$score, b$score, c$score, d$score, e$score,
        col = "#4D648D", xlab="Number of Exon", ylab="Splicing Efficiency", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), main="")
boxplot(a1$score, b1$score, c1$score, d1$score, e1$score,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()