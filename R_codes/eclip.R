#NMD analysis and plots with eclip data

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

data_pcg <- read.csv("mrna_eclip_h.csv", header = FALSE)
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

#scatterplot
tiff("/home/dell/Documents/koushiki/final_fig/eclip_fig_1.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5))
plot(data_pcg$Transcript, data_pcg$Exon, xlim=c(0,90), ylim = c(0,400),
     xlab = "Number of transcripts", ylab = "Number of exons",
     main="",col="#E8A735", pch=20, font.axis=2, cex.axis=1, font.lab=2, cex.lab=1.2 )
abline(lm(data_pcg$Exon~data_pcg$Transcript), lty="solid", col="#EB8A3E", lwd=3)
points(data_lncRNA$Transcript, data_lncRNA$Exon, col="#4D648D",pch=20)
abline(lm(data_lncRNA$Exon~data_lncRNA$Transcript), lty="solid", col="midnightblue", lwd=3)
dev.off()

#boxplot
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
tiff("/home/dell/Documents/koushiki/final_fig/eclip_fig_2.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(a$Transcript, b$Transcript, c$Transcript, d$Transcript, e$Transcript,
        col = "#4D648D", xlab="Number of Exon", ylab="Number of Transcript", boxwex=0.9,xaxs="i",xaxt='n',
        at=c(1,4,7,10,13), outline=F, xlim=c(0,15), ylim=c(0,50), main="")
boxplot(a1$Transcript, b1$Transcript, c1$Transcript, d1$Transcript, e1$Transcript,
        add= TRUE,col = "#E8A735", boxwex=0.9, at=c(2,5,8,11,14), outline=F,axes=F)
axis(1, at=c(1.5,4.5,7.5,10.5,13.5), labels=c("0-10","11-20","21-50","51-100",">100"),font=2,cex=0.9)
dev.off()


#analysis based on mean TPE and plots
data_lncRNA <- read.csv("lncrna_main_h.csv", header = FALSE)
data_lncRNA <- na.omit(data_lncRNA)
names(data_lncRNA)
colnames(data_lncRNA) <- c("Gene","Gene_length","Transcript","Exon")
summary(data_lncRNA)
str(data_lncRNA)
data_lncRNA$Gene_length <- as.numeric(data_lncRNA$Gene_length)
data_lncRNA$Transcript <- as.numeric(data_lncRNA$Transcript)
data_lncRNA$Exon <- as.numeric(data_lncRNA$Exon)
data_lncRNA <- mutate(data_lncRNA, TPE=Transcript/Exon)
median(data_lncRNA$TPE)

data_pcg <- read.csv("mrna_eclip_h.csv", header = FALSE)
data_pcg <- na.omit(data_pcg)
names(data_pcg)
colnames(data_pcg) <- c("Gene","Gene_length","Transcript","Exon")
summary(data_pcg)
str(data_pcg)
data_pcg$Gene_length <- as.numeric(data_pcg$Gene_length)
data_pcg$Transcript <- as.numeric(data_pcg$Transcript)
data_pcg$Exon <- as.numeric(data_pcg$Exon)
data_pcg <- mutate(data_pcg, TPE=Transcript/Exon)
median(data_pcg$TPE)

#boxplot
tiff("/home/dell/Documents/koushiki/final_fig/eclip_fig_3.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(0.2,5,0.2,0.2))
boxplot(data_lncRNA$TPE, data_pcg$TPE, col = c("#4D648D", "#E8A735"),ylab="TPE", boxwex=2.0,
        xaxs="i",at=c(1,4), cex.axis=1.5, font.axix=2,cex.lab=1.8, font.lab=2,
        outline=F, xlim=c(0,5), ylim=c(0.0,0.8), main="")
dev.off()

#density plot
data_lncRNA1 <- read.csv("lncrna_main_h.csv", header = F)
colnames(data_lncRNA1) <- c('Gene','Gene_length','Transcript','Exon')
data_CDS1 <- read.csv("mrna_eclip_h.csv", header = F)
colnames(data_CDS1) <- c('Gene','Gene_length','Transcript','Exon')

library(dplyr)
data_lncRNA_new1 <- mutate(data_lncRNA1, TPE = Transcript/Exon)
hist(data_lncRNA_new1$TPE)
data_CDS_new1 <- mutate(data_CDS1, TPE = Transcript/Exon)
hist(data_CDS_new1$TPE)
mean(data_lncRNA_new1$TPE)
mean(data_CDS_new1$TPE, na.rm=T)
sd(data_lncRNA_new1$TPE)
sd(data_CDS_new1$TPE, na.rm=T)
length(data_lncRNA_new1$TPE)
length(data_CDS_new1$TPE)

h3 <- density(data_CDS_new1$TPE, na.rm = T, kernel = "cosine", bw = "bcv", adjust = 3)
h4 <- density(data_lncRNA_new1$TPE, kernel = "cosine", bw = "bcv", adjust = 3)
tiff("/home/dell/Documents/koushiki/final_fig/eclip_fig_4.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4.1,4.5,1.5,1.5))
plot(h3, xlim=c(0,1), col="#E8A735", ylim=c(0,3),xlab="Transcript Per Exon", lwd = 4, main="",cex.lab = 1.3, cex.axis=0.8, font.lab=2, font.axis=2)
lines(h4, col="#4D648D", lwd=4)
dev.off()