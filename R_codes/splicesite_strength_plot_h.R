#plot of splicesite strength

# 5'ss strength (MAXENT, MDD, MM, WMM)
data_lncRNA_low <- read.table("final_lncrna_low_5sss_h.txt", sep = "\t")
colnames(data_lncRNA_low) <- c("motif","MAXENT","MDD","MM","WMM")
head(data_lncRNA_low)
data_lncRNA_low$motif <- as.character(data_lncRNA_low$motif)
summary(data_lncRNA_low)

data_lncRNA_high <- read.table("final_lncrna_high_5sss_h.txt", sep = "\t")
colnames(data_lncRNA_high) <- c("motif","MAXENT","MDD","MM","WMM")
head(data_lncRNA_high)
data_lncRNA_high$motif <- as.character(data_lncRNA_high$motif)
summary(data_lncRNA_high)

data_CDS_low <- read.table("final_mrna_low_5sss_h.txt", sep = "\t")
colnames(data_CDS_low) <- c("motif","MAXENT","MDD","MM","WMM")
head(data_CDS_low)
data_CDS_low$motif <- as.character(data_CDS_low$motif)
summary(data_CDS_low)

data_CDS_high <- read.table("final_mrna_high_5sss_h.txt", sep = "\t")
colnames(data_CDS_high) <- c("motif","MAXENT","MDD","MM","WMM")
head(data_CDS_high)
data_CDS_high$motif <- as.character(data_CDS_high$motif)
summary(data_CDS_high)

tiff("/home/dell/Documents/koushiki/final_fig/maxentscan_5_h.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA_low$MAXENT, data_lncRNA_low$MDD, data_lncRNA_low$MM, data_lncRNA_low$WMM,
        col = "#3A5199", ylab="5' Splice site strength", xlab="Scoring models", boxwex=0.8,xaxs="i",xaxt='n',
        at=c(1,6,11,16), outline=F, xlim=c(0,20), ylim=c(0,20),main="")
boxplot(data_lncRNA_high$MAXENT, data_lncRNA_high$MDD, data_lncRNA_high$MM, data_lncRNA_high$WMM,
        add= TRUE,col = "#31A2AC", boxwex=0.8, at=c(2,7,12,17), outline=F,axes=F)
boxplot(data_CDS_low$MAXENT, data_CDS_low$MDD, data_CDS_low$MM, data_CDS_low$WMM,
        add= TRUE,col = "#CB0000", boxwex=0.8, at=c(3,8,13,18), outline=F,axes=F)
boxplot(data_CDS_high$MAXENT, data_CDS_high$MDD, data_CDS_high$MM, data_CDS_high$WMM,
        add= TRUE,col = "#EB8A3E", boxwex=0.8, at=c(4,9,14,19), outline=F,axes=F)
axis(1, at=c(2.5,7.5,12.5,17.5), labels=c("MAXENT","MDD","MM","WMM"),font=2,cex=0.9)
dev.off()

#ks test
#MAXENT
ks.test(data_lncRNA_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_low$MAXENT, data_CDS_low$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_low$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_high$MAXENT, data_CDS_low$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_high$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
ks.test(data_CDS_low$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
#MDD
ks.test(data_lncRNA_low$MDD, data_lncRNA_high$MDD, alternative = "two.sided")
ks.test(data_lncRNA_low$MDD, data_CDS_low$MDD, alternative = "two.sided")
ks.test(data_lncRNA_low$MDD, data_CDS_high$MDD, alternative = "two.sided")
ks.test(data_lncRNA_high$MDD, data_CDS_low$MDD, alternative = "two.sided")
ks.test(data_lncRNA_high$MDD, data_CDS_high$MDD, alternative = "two.sided")
ks.test(data_CDS_low$MDD, data_CDS_high$MDD, alternative = "two.sided")
#MM
ks.test(data_lncRNA_low$MM, data_lncRNA_high$MM, alternative = "two.sided")
ks.test(data_lncRNA_low$MM, data_CDS_low$MM, alternative = "two.sided")
ks.test(data_lncRNA_low$MM, data_CDS_high$MM, alternative = "two.sided")
ks.test(data_lncRNA_high$MM, data_CDS_low$MM, alternative = "two.sided")
ks.test(data_lncRNA_high$MM, data_CDS_high$MM, alternative = "two.sided")
ks.test(data_CDS_low$MM, data_CDS_high$MM, alternative = "two.sided")
#WMM
ks.test(data_lncRNA_low$WMM, data_lncRNA_high$WMM, alternative = "two.sided")
ks.test(data_lncRNA_low$WMM, data_CDS_low$WMM, alternative = "two.sided")
ks.test(data_lncRNA_low$WMM, data_CDS_high$WMM, alternative = "two.sided")
ks.test(data_lncRNA_high$WMM, data_CDS_low$WMM, alternative = "two.sided")
ks.test(data_lncRNA_high$WMM, data_CDS_high$WMM, alternative = "two.sided")
ks.test(data_CDS_low$WMM, data_CDS_high$WMM, alternative = "two.sided")

#wilcox test
#MAXENT
wilcox.test(data_lncRNA_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_lncRNA_low$MAXENT, alternative = "g")
wilcox.test(data_CDS_high$MAXENT, data_lncRNA_low$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_high$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_CDS_high$MAXENT, alternative = "g")
#MDD
wilcox.test(data_lncRNA_low$MDD, data_lncRNA_high$MDD, alternative = "g")
wilcox.test(data_CDS_low$MDD, data_lncRNA_low$MDD, alternative = "g")
wilcox.test(data_CDS_high$MDD, data_lncRNA_low$MDD, alternative = "g")
wilcox.test(data_CDS_low$MDD, data_lncRNA_high$MDD, alternative = "g")
wilcox.test(data_CDS_high$MDD, data_lncRNA_high$MDD, alternative = "g")
wilcox.test(data_CDS_low$MDD, data_CDS_high$MDD, alternative = "g")
#MM
wilcox.test(data_lncRNA_low$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_low$MM, data_lncRNA_low$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_lncRNA_low$MM, alternative = "g")
wilcox.test(data_CDS_low$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_CDS_low$MM, alternative = "g")
#WMM
wilcox.test(data_lncRNA_low$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_low$WMM, data_lncRNA_low$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_lncRNA_low$WMM, alternative = "g")
wilcox.test(data_CDS_low$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_CDS_low$WMM, alternative = "g")

# 3'ss strength (MAXENT, MM, WMM)
data_lncRNA_low <- read.table("final_lncrna_low_3sss_h.txt", sep = "\t")
colnames(data_lncRNA_low) <- c("motif","MAXENT","MM","WMM")
head(data_lncRNA_low)
data_lncRNA_low$motif <- as.character(data_lncRNA_low$motif)
summary(data_lncRNA_low)

data_lncRNA_high <- read.table("final_lncrna_high_3sss_h.txt", sep = "\t")
colnames(data_lncRNA_high) <- c("motif","MAXENT","MM","WMM")
head(data_lncRNA_high)
data_lncRNA_high$motif <- as.character(data_lncRNA_high$motif)
summary(data_lncRNA_high)

data_CDS_low <- read.table("final_mrna_low_3sss_h.txt", sep = "\t")
colnames(data_CDS_low) <- c("motif","MAXENT","MM","WMM")
head(data_CDS_low)
data_CDS_low$motif <- as.character(data_CDS_low$motif)
summary(data_CDS_low)

data_CDS_high <- read.table("final_mrna_high_3sss_h.txt", sep = "\t")
colnames(data_CDS_high) <- c("motif","MAXENT","MM","WMM")
head(data_CDS_high)
data_CDS_high$motif <- as.character(data_CDS_high$motif)
summary(data_CDS_high)

tiff("/home/dell/Documents/koushiki/final_fig/maxentscan_3_h.tiff",width = 1200, height = 1200, res = 300)
par(mar=c(4,4,0.5,0.5),cex.axis=0.9, font.axis=2,cex.lab=1.3, font.lab=2)
boxplot(data_lncRNA_low$MAXENT, data_lncRNA_low$MM, data_lncRNA_low$WMM,
        col = "#3A5199", ylab="3' Splice site strength", xlab="Scoring models", boxwex=0.8,xaxs="i",xaxt='n',
        at=c(1,6,11), outline=F, xlim=c(0,15), ylim=c(0,20),main="")
boxplot(data_lncRNA_high$MAXENT, data_lncRNA_high$MM, data_lncRNA_high$WMM,
        add= TRUE,col = "#31A2AC", boxwex=0.8, at=c(2,7,12), outline=F,axes=F)
boxplot(data_CDS_low$MAXENT, data_CDS_low$MM, data_CDS_low$WMM,
        add= TRUE,col = "#CB0000", boxwex=0.8, at=c(3,8,13), outline=F,axes=F)
boxplot(data_CDS_high$MAXENT, data_CDS_high$MM, data_CDS_high$WMM,
        add= TRUE,col = "#EB8A3E", boxwex=0.8, at=c(4,9,14), outline=F,axes=F)
axis(1, at=c(2.5,7.5,12.5), labels=c("MAXENT","MM","WMM"),font=2,cex=0.9)
dev.off()

#ks test
#MAXENT
ks.test(data_lncRNA_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_low$MAXENT, data_CDS_low$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_low$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_high$MAXENT, data_CDS_low$MAXENT, alternative = "two.sided")
ks.test(data_lncRNA_high$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
ks.test(data_CDS_low$MAXENT, data_CDS_high$MAXENT, alternative = "two.sided")
#MM
ks.test(data_lncRNA_low$MM, data_lncRNA_high$MM, alternative = "two.sided")
ks.test(data_lncRNA_low$MM, data_CDS_low$MM, alternative = "two.sided")
ks.test(data_lncRNA_low$MM, data_CDS_high$MM, alternative = "two.sided")
ks.test(data_lncRNA_high$MM, data_CDS_low$MM, alternative = "two.sided")
ks.test(data_lncRNA_high$MM, data_CDS_high$MM, alternative = "two.sided")
ks.test(data_CDS_low$MM, data_CDS_high$MM, alternative = "two.sided")
#WMM
ks.test(data_lncRNA_low$WMM, data_lncRNA_high$WMM, alternative = "two.sided")
ks.test(data_lncRNA_low$WMM, data_CDS_low$WMM, alternative = "two.sided")
ks.test(data_lncRNA_low$WMM, data_CDS_high$WMM, alternative = "two.sided")
ks.test(data_lncRNA_high$WMM, data_CDS_low$WMM, alternative = "two.sided")
ks.test(data_lncRNA_high$WMM, data_CDS_high$WMM, alternative = "two.sided")
ks.test(data_CDS_low$WMM, data_CDS_high$WMM, alternative = "two.sided")



#wilcox test
#MAXENT
wilcox.test(data_lncRNA_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_lncRNA_low$MAXENT, alternative = "g")
wilcox.test(data_CDS_high$MAXENT, data_lncRNA_low$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_high$MAXENT, data_lncRNA_high$MAXENT, alternative = "g")
wilcox.test(data_CDS_low$MAXENT, data_CDS_high$MAXENT, alternative = "g")
#MM
wilcox.test(data_lncRNA_low$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_low$MM, data_lncRNA_low$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_lncRNA_low$MM, alternative = "g")
wilcox.test(data_CDS_low$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_lncRNA_high$MM, alternative = "g")
wilcox.test(data_CDS_high$MM, data_CDS_low$MM, alternative = "g")
#WMM
wilcox.test(data_lncRNA_low$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_low$WMM, data_lncRNA_low$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_lncRNA_low$WMM, alternative = "g")
wilcox.test(data_CDS_low$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_lncRNA_high$WMM, alternative = "g")
wilcox.test(data_CDS_high$WMM, data_CDS_low$WMM, alternative = "g")


