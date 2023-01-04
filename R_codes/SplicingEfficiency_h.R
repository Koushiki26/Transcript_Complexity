#splicing efficiency
data_lncRNA <- read.csv("lncrna_main_h.csv", header = F)
head(data_lncRNA)
colnames(data_lncRNA) <- c("gene_ID", "gene_length", "transcript","exon")
data_lncRNA$TC <- data_lncRNA$transcript/data_lncRNA$exon
data_CDS <- read.csv("mrna_main_h.csv", header = F)
head(data_CDS)
colnames(data_CDS) <- c("gene_ID", "gene_length", "transcript","exon")
data_CDS$TC <- data_CDS$transcript/data_CDS$exon
K562_data <- read.table("K562_output.tsv", header = T, sep="\t")
head(K562_data)
summary(K562_data)
A549_data <- read.table("A549_output.tsv", header = T, sep = "\t")
head(A549_data)
summary(A549_data)

library(dplyr)
K562_data$TC <- ''
K562_data$exon <- ''
for(i in 1:length(K562_data$score))
{
  TC <- data_lncRNA %>% filter(as.character(gene_ID) == as.character(K562_data$gene_ID[i]))
  K562_data$TC[i] = TC$TC[1]
  K562_data$exon[i] = TC$exon[1]
}
K562_data <- na.omit(K562_data)
write.csv(K562_data, file = "lncRNA_K562.csv", col.names = T, row.names = F)

K562_data$TC <- ''
K562_data$exon <- ''
for(i in 1:length(K562_data$score))
{
  TC <- data_CDS %>% filter(as.character(gene_ID) == as.character(K562_data$gene_ID[i]))
  K562_data$TC[i] = TC$TC[1]
  K562_data$exon[i] = TC$exon[1]
}
K562_data <- na.omit(K562_data)
write.csv(K562_data, file = "mRNA_K562.csv", col.names = T, row.names = F)

data_lncRNA <- read.csv("lncRNA_K562.csv", header = T)
head(data_lncRNA)
data_pcg <- read.csv("mRNA_K562.csv", header = T)
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
tiff("/home/dell/Documents/koushiki/final_fig/SE_bin1_K562.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$score, b$score, c$score, d$score, e$score,
        col = "#4D648D", xlab="Number of Exon", ylab="Splicing Efficiency", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), main="")
boxplot(a1$score, b1$score, c1$score, d1$score, e1$score,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()

A549_data$TC <- ''
A549_data$exon <- ''
for(i in 1:length(A549_data$score))
{
  TC <- data_lncRNA %>% filter(as.character(gene_ID) == as.character(A549_data$gene_ID[i]))
  A549_data$TC[i] = TC$TC[1]
  A549_data$exon[i] = TC$exon[1]
}
A549_data <- na.omit(A549_data)
write.csv(A549_data, file = "lncRNA_A549.csv", col.names = T, row.names = F)

A549_data$TC <- ''
A549_data$exon <- ''
for(i in 1:length(A549_data$score))
{
  TC <- data_CDS %>% filter(as.character(gene_ID) == as.character(A549_data$gene_ID[i]))
  A549_data$TC[i] = TC$TC[1]
  A549_data$exon[i] = TC$exon[1]
}
A549_data <- na.omit(A549_data)
write.csv(A549_data, file = "mRNA_A549.csv", col.names = T, row.names = F)

data_lncRNA <- read.csv("lncRNA_A549.csv", header = T)
head(data_lncRNA)
data_pcg <- read.csv("mRNA_A549.csv", header = T)
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
tiff("/home/dell/Documents/koushiki/final_fig/SE_bin1_A549.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$score, b$score, c$score, d$score, e$score,
        col = "#4D648D", xlab="Number of Exon", ylab="Splicing Efficiency", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), main="")
boxplot(a1$score, b1$score, c1$score, d1$score, e1$score,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()

