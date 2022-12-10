#correlation between the number of transcript and number of exon

data_lncRNA <- read.csv("lncrna_main_h.csv", header = FALSE)
head(data_lncRNA)
names(data_lncRNA)
colnames(data_lncRNA)[1] <- "Gene"
colnames(data_lncRNA)[2] <- "Gene_length"
colnames(data_lncRNA)[3] <- "Transcript"
colnames(data_lncRNA)[4] <- "Exon"
summary(data_lncRNA)
str(data_lncRNA)
data_lncRNA$Gene_length <- as.numeric(data_lncRNA$Gene_length)
data_lncRNA$Transcript <- as.numeric(data_lncRNA$Transcript)
data_lncRNA$Exon <- as.numeric(data_lncRNA$Exon)
length(data_lncRNA$Transcript)
length(data_lncRNA$Exon)
hist(data_lncRNA$Transcript)
hist(data_lncRNA$Exon)
cor.test(data_lncRNA$Transcript, data_lncRNA$Exon, method = c('spearman'))
cor.test(data_lncRNA$Exon, data_lncRNA$Gene_length, method = c('spearman'))
cor.test(data_lncRNA$Transcript, data_lncRNA$Gene_length, method = c('spearman'))
cor.test(data_lncRNA$Transcript, data_lncRNA$Exon, method = c('pearson'))

data_pcg <- read.csv("mrna_main_h.csv", header = FALSE)
head(data_pcg)
names(data_pcg)
colnames(data_pcg)[1] <- "Gene"
colnames(data_pcg)[2] <- "Gene_length"
colnames(data_pcg)[3] <- "Transcript"
colnames(data_pcg)[4] <- "Exon"
summary(data_pcg)
str(data_pcg)
data_pcg$Gene_length <- as.numeric(data_pcg$Gene_length)
data_pcg$Transcript <- as.numeric(data_pcg$Transcript)
data_pcg$Exon <- as.numeric(data_pcg$Exon)
length(data_pcg$Transcript)
length(data_pcg$Exon)
hist(data_pcg$Transcript)
hist(data_pcg$Exon)
cor.test(data_pcg$Transcript, data_pcg$Exon, method = c('spearman'))
cor.test(data_pcg$Exon, data_pcg$Gene_length, method = c('spearman'))
cor.test(data_pcg$Transcript, data_pcg$Gene_length, method = c('spearman'))
cor.test(data_pcg$Transcript, data_pcg$Exon, method = c('pearson'))

#scatterplot
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_1a1.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5))
plot(data_pcg$Transcript, data_pcg$Exon, xlim=c(0,90), ylim = c(0,400),
     xlab = "Number of transcripts", ylab = "Number of exons",
     main="",col="#E8A735", pch=20, font.axis=2, cex.axis=1, font.lab=2, cex.lab=1.2 )
abline(lm(data_pcg$Exon~data_pcg$Transcript), lty="solid", col="#EB8A3E", lwd=3)
points(data_lncRNA$Transcript, data_lncRNA$Exon, col="#4D648D",pch=20)
abline(lm(data_lncRNA$Exon~data_lncRNA$Transcript), lty="solid", col="midnightblue", lwd=3)
dev.off()

#boxplot
library(dplyr)
a <- data_lncRNA %>% filter(Exon<11)
count(a)
b <- data_lncRNA %>% filter(Exon>10&Exon<21)
count(b)
c <- data_lncRNA %>% filter(Exon>20&Exon<51)
count(c)
d <- data_lncRNA %>% filter(Exon>50&Exon<101)
count(d)
e <- data_lncRNA %>% filter(Exon>100)
count(e)
a1 <- data_pcg %>% filter(Exon<11)
count(a1)
b1 <- data_pcg %>% filter(Exon>10&Exon<21)
count(b1)
c1 <- data_pcg %>% filter(Exon>20&Exon<51)
count(c1)
d1 <- data_pcg %>% filter(Exon>50&Exon<101)
count(d1)
e1 <- data_pcg %>% filter(Exon>100)
count(e1)
tiff("/home/dell/Documents/koushiki/final_fig/supplementary_fig_1b1.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$Transcript, b$Transcript, c$Transcript, d$Transcript, e$Transcript,
        col = "#4D648D", xlab="Number of Exon", ylab="Number of Transcript", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), ylim=c(0,50), main="")
boxplot(a1$Transcript, b1$Transcript, c1$Transcript, d1$Transcript, e1$Transcript,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()
